/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 **/

#ifndef COMMON_PROGRAM_HPP
#define COMMON_PROGRAM_HPP

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>

#include "common/Debug.hpp"
#include "common/Exceptions.hpp"
#include "common/Version.hpp"

namespace dragenos {
namespace common {

namespace bpo = boost::program_options;

/**
 ** Encapsulation of the pocessing of the command line options.
 **
 ** TODO: add config file and environment options
 **/
class Options : boost::noncopyable {
public:
  enum Action { RUN, HELP, VERSION, ABORT };
  Options();
  virtual ~Options() {}
  Action             parse(int argc, const char* const argv[]);
  std::string        usage() const;
  int                argc() const { return argc_; }
  const char* const* argv() const { return argv_; }

  std::string getCommandLine() const
  {
    const std::vector<std::string> allOptions(argv(), argv() + argc());
    return boost::join(allOptions, " ");
  }

  bool exists(const std::string& opt) const { return vm_.count(opt); }

protected:
  bool                                version_ = false;
  bpo::options_description            namedOptions_;
  bpo::options_description            unnamedOptions_;
  bpo::positional_options_description positionalOptions_;

  typedef boost::shared_ptr<boost::program_options::option_description> OptionDescriptionPtr;
  typedef std::vector<OptionDescriptionPtr>                             OptionDescriptionPtrs;
  std::string helpDefaults(const OptionDescriptionPtrs& options) const;
  std::string help(const OptionDescriptionPtrs& options, const bool markdown) const;

private:
  virtual std::string usagePrefix() const = 0;
  virtual std::string usageSuffix() const { return ""; }
  virtual void        postProcess(bpo::variables_map&) {}

  static const unsigned MARKDOWN_LINE_LENGTH = 120;
  // "parse" will store the state in vm_ so that "usage" can access the details of parsed command line
  bpo::variables_map vm_;
  int                argc_ = 0;
  const char* const* argv_ = 0;
};

/**
 ** Unified behavior of all programs.
 **/
template <class O>
void run(void (*callback)(const O&), int argc, char* argv[])
{
  // fix for Failed to parse the options: locale::facet::_S_create_c_locale name not valid
  setenv("LC_ALL", "C", 0);
  try {
    O                        options;
    const typename O::Action action = options.parse(argc, argv);
    if (O::RUN == action) {
      callback(options);
    } else if (O::HELP == action) {
      std::cout << options.usage() << std::endl;
      exit(0);
    } else if (O::VERSION == action) {
      std::cout << Version::string() << std::endl;
      exit(0);
    } else {
      //            std::clog << options.usage() << std::endl;
      exit(1);
    }
  } catch (const dragenos::common::ExceptionData& exception) {
    std::clog << "Error: " << exception.getContext() << ": " << exception.getMessage() << std::endl;
    exit(1);
  } catch (const boost::exception& e) {
    std::clog << "Error: boost::exception: " << boost::diagnostic_information(e) << std::endl;
    exit(2);
  } catch (const std::exception& e) {
    std::clog << e.what() << std::endl;
    exit(3);
  }
}

}  // namespace common
}  // namespace dragenos

#endif  // #ifndef COMMON_PROGRAM_HPP
