#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

#set -o xtrace

this_file="$( readlink -f "${BASH_SOURCE[0]}" )"
this_file_name=$( basename $this_file )
meta_dir="$( dirname $this_file )"
repo_dir="$( dirname $meta_dir )"
src_dir=$repo_dir/src
source $meta_dir/common.bash
n_cpu=$(nproc)
n_jobs=$(( n_cpu > 16 ? 16 : n_cpu ))
show_help(){
    echo "$this_file_name [OPTIONS] [FILE...]"
    echo "-h --help          show this help menu"
    echo "-c --check         only check for formatting changes without actually formatting"
    echo "-s --selection     selection of files for formatting (only files ending in .cpp .cc .c .hpp .hh and .h will be processed unless --unfiltered flag is specified). possible values:"
    echo "                     -        read/write source code from/to standard input/standard output"
    echo "                     staged   staged files"
    echo "                     changed  staged and changed files relative to HEAD (includes untracked files)"
    echo "                     master   staged and changed files relative to the most recent shared commit with origin/master (default)"
    echo "                     all      files in this repository src directory"
    echo "                     files    use list of files from standard input (must be NUL separated)"
    echo "-u --unfiltered    reformat all specified files regardless of any filename extension"
    echo "-n --num-jobs      number of files to process in parallel (default $n_jobs)"
}

OPTS=`getopt -o hcs:n: --long help,check,selection:,num-jobs: -n $this_file_name -- "$@"`
eval set -- "$OPTS"
check=false
selection=args
unfiltered=false
while true; do
  case "$1" in
    -h | --help      )  show_help; exit 0 ;;
    -c | --check     )  check=true ; shift ;;
    -s | --selection )
       case "$2" in
         -       ) unfiltered=true ;&
         staged  ) ;&
         changed ) ;&
         master  ) ;&
         all     ) ;&
         files   )
                        selection=$2 ; shift 2 ;;
         *       )      echo "Error: invalid value for option $1: $2" ; echo ; show_help ; exit 1 ;;
       esac ;;
    -u | --unfiltered ) unfiltered=true ; shift ;;
    -n | --num-jobs  )
                        re='^[1-9]+[0-9]*$'
                        if ! [[ "$2" =~ $re ]] ; then
                          echo "Error: value for option '$1' must be an integer greater than 0" >&2; exit 1
                        fi
                        n_jobs=$2 ; shift 2 ;;
     --              )  shift ; break ;;
     *               )  echo "Internal error" ; exit 1 ;;
  esac
done
files="$@"
if [ "$selection" != "args" ] && [ "$files" != "" ] ; then
  echo "Error: selection option '$selection' cannot be combined with FILE arguments"
  exit 1
fi
if [ "$selection" == "args" ] && [ "$files" == "" ] ; then
  selection=master
fi

get_file_args(){
  for f in $files;
  do
  printf '%s\0' "$f"
  done
}

get_repo_files() {
  ( cd $src_dir && find_files_current_dir ; )
}

clang_format=$meta_dir/clang-format
if [[ ${OS:=} == Windows* ]] ; then
  clang_format=$meta_dir/clang-format.exe
fi

get_files() {
  case "$1" in
    staged  ) cd $repo_dir ; git_staged_files ;;
    changed ) cd $repo_dir ; get_files staged ; git_unstaged_and_untracked_files ;;
    master  ) cd $repo_dir ; get_files changed ; git_committed_since_master_files ;;
    all     ) get_repo_files ;;
    -       ) ;&
    files   ) cat ;;
    args    ) get_file_args ;;
  esac
}

get_source_files() {
  if [ "$unfiltered" = true ] ; then
    get_files $selection
  else
    get_files $selection | get_c_and_cpp_src_and_headers
  fi
}

print_path_to_stderr(){
  echo $1 > /dev/stderr
}

check_changes() {
  if [ "$check" = true ] ; then
    grep '<replacement ' > /dev/null && exitCode=$? || exitCode=$?
    if [ "$exitCode" -eq "0" ] ; then
      print_path_to_stderr $1
    fi
  else
    cat
  fi
}

run_clang_format() {
 if [ "$check" = true ] ; then
   if [ "$selection" = "-" ] ; then
     $clang_format -style=file -assume-filename=$1 -output-replacements-xml
   else
     $clang_format -style=file -output-replacements-xml $1
   fi
 elif [ "$selection" = "-" ] ; then
    $clang_format -style=file -assume-filename=$1
 else
    $clang_format -style=file -i $1
 fi
}

fix_line_endings() {
  if [ "$check" = true ] ; then
    if [ "$selection" = "-" ] ; then
      sed -b '/\r$/{q1}' && exitCode=$? || exitCode=$?
    else
      sed -b '/\r$/{q1}' "$1" > /dev/null && exitCode=$? || exitCode=$?
    fi
    if [[ "$exitCode" != "0" ]] ; then
      print_path_to_stderr $1
    fi
  elif [ "$selection" = "-" ] ; then
    sed -b 's/\r$//'
  else
    sed -bi 's/\r$//' $1
  fi
}

reformat_and_check() {
  set -o errexit
  set -o nounset
  set -o pipefail
  case "$selection" in
    staged  ) cd $repo_dir ;;
    changed ) cd $repo_dir ;;
    master  ) cd $repo_dir ;;
    all     ) cd $src_dir ;;
  esac
  fix_line_endings $1 | run_clang_format $1 | check_changes $1
}

export -f reformat_and_check
export -f fix_line_endings
export -f check_changes
export -f run_clang_format
export -f print_path_to_stderr
export selection
export check
export clang_format
export src_dir
export repo_dir

if [ "$selection" = "-" ] ; then
    (reformat_and_check $src_dir/dummy.cpp)
fi

set +o errexit
files_need_changing=$( (get_source_files | xargs --no-run-if-empty -0 -P$n_jobs -I % bash -c "reformat_and_check %") 2>&1)
files_need_changing=$(echo $files_need_changing | xargs -n1 | sort --unique)
if [[ "$?" != "0" ]] ; then
    echo "Internal Error!"
    echo $files_need_changing
    exit 1;
fi
set -o errexit

if [[ "$check" != true ]]; then exit 0 ; fi

if [[ "$files_need_changing" == "" ]] ; then
  echo "Formatting looks good!"
  exit 0
else
  echo "Code needs formatting!"
  if [[ "$selection" != "-" ]] ; then
    echo "Files not matching expected format:"
    for f in $files_need_changing ;
    do
    if [[ "$selection" == "all" ]] ; then
      f=$(echo $f | sed -e 's/^\./src/')
    fi
    echo "$f"
    done
  fi
  exit 1
fi
