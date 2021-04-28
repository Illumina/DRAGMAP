exclude_dirs_find="-type d ( -path ./third_party -o -path ./workbench -o -path ./host/cram -o -path ./host/repeat_genotyping/expansion_hunter ) -prune"
exclude_dirs_grep='src/third_party|src/workbench|src/host/cram|src/host/repeat_genotyping/expansion_hunter'
exclude_dirs_grep='${exclude_dirs_grep}|src/host/dipc/.*pb.cc|src/host/dipc/.*pb.h'

find_files_current_dir() {
 \find . $exclude_dirs_find -o -type f -print0
}

exclude_dirs() {
 \grep -zvE $exclude_dirs_grep || true
}

git_staged_files() {
 git --no-pager diff -z --cached --name-only --diff-filter=ACMRT | exclude_dirs
}

git_unstaged_and_untracked_files() {
 git --no-pager ls-files -z -m -o --exclude-standard | exclude_dirs
}

git_committed_since_master_files() {
 git --no-pager diff -z --name-only --diff-filter=ACMRT origin/master... | exclude_dirs
}

get_c_and_cpp_src_and_headers() {
  \grep --null-data -E '(\.cpp|\.cc|\.c|\.hpp|\.hh|\.h)$' || true
}
