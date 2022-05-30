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

#include <fstream>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "common/Program.hpp"

namespace dragenos {
namespace common {

#define HELP_STR "help"
#define HELP_DEFAULTS_STR "help-defaults"

Options::Options()
{
  namedOptions_.add_options()(HELP_STR ",h", "produce help message and exit");
  namedOptions_.add_options()(
      HELP_DEFAULTS_STR, "produce tab-delimited list of command line options and their default values");
  namedOptions_.add_options()("version,V", bpo::bool_switch(&version_), "print program version information");
  namedOptions_.add_options()(
      "response-file", bpo::value<std::string>(), "file with more command line arguments");
}

Options::Action Options::parse(int argc, const char* const argv[])
{
  argc_ = argc;
  argv_ = argv;
  try {
    bpo::options_description allOptions("Allowed options");
    allOptions.add(namedOptions_).add(unnamedOptions_);
    vm_.clear();
    bpo::store(
        bpo::command_line_parser(argc, argv).options(allOptions).positional(positionalOptions_).run(), vm_);

    if (vm_.count("response-file")) {
      // Load the file and tokenize it
      std::ifstream ifs(vm_["response-file"].as<std::string>().c_str());
      if (!ifs) {
        std::clog << "Could not open response file: " << vm_["response-file"].as<std::string>().c_str()
                  << std::endl;
        return ABORT;
      }
      // Read the whole file into a string
      std::stringstream ss;
      ss << ifs.rdbuf();
      if (!ifs) {
        std::clog << "Could not read response file: " << vm_["response-file"].as<std::string>().c_str()
                  << std::endl;
        return ABORT;
      }
      // Split the file content
      boost::char_separator<char>                   sep(" \n\r");
      std::string                                   ResponsefileContents(ss.str());
      boost::tokenizer<boost::char_separator<char>> tok(ResponsefileContents, sep);
      std::vector<std::string>                      args;
      copy(tok.begin(), tok.end(), std::back_inserter(args));
      // Parse the file and store the options
      bpo::store(bpo::command_line_parser(args).options(allOptions).run(), vm_);
    }

    bpo::notify(vm_);
    if (vm_.count(HELP_STR) || vm_.count(HELP_DEFAULTS_STR)) {
      return HELP;
    } else if (version_) {
      return VERSION;
    }
    postProcess(vm_);
  } catch (const boost::program_options::multiple_values& e) {
    std::clog << usage() << std::endl;
    std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
    return ABORT;
  } catch (const boost::program_options::multiple_occurrences& e) {
    std::clog << usage() << std::endl;
    std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
    return ABORT;
  } catch (const boost::program_options::required_option& e) {
    std::clog << usage() << std::endl;
    std::clog << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
    return ABORT;
  } catch (const boost::exception& e) {
    std::clog << usage() << std::endl;
    std::clog << "Failed to parse the options: " << boost::diagnostic_information(e) << std::endl;
    return ABORT;
  } catch (const std::exception& e) {
    std::clog << usage() << std::endl;
    std::clog << "Failed to parse the options: " << e.what() << std::endl;
    return ABORT;
  }
  return RUN;
}

bool compareOptionName(
    const boost::shared_ptr<boost::program_options::option_description>& left,
    const boost::shared_ptr<boost::program_options::option_description>& right)
{
  return left->long_name() < right->long_name();
}

std::string Options::helpDefaults(const OptionDescriptionPtrs& sortedOptions) const
{
  std::string ret;
  BOOST_FOREACH (const OptionDescriptionPtr& odPtr, sortedOptions) {
    ret += odPtr->long_name() + "\t" + odPtr->format_parameter() + "\n";
  }
  return ret;
}

static unsigned short getTerminalColumns()
{
  unsigned short int ws_row = 0;
  unsigned short int ws_col = 0;

  if (-1 == getTerminalWindowSize(ws_row, ws_col)) {
    return bpo::options_description::m_default_line_length;
  }

  return ws_col;
}

std::string Options::help(const OptionDescriptionPtrs& sortedOptions) const
{
  std::ostringstream       os;
  const unsigned           effectiveLineLength = getTerminalColumns();
  bpo::options_description printedDescriptions(
      "Command line options",
      effectiveLineLength < 50 ? bpo::options_description::m_default_line_length : effectiveLineLength,
      effectiveLineLength < 50 ? bpo::options_description::m_default_line_length / 2
                               : effectiveLineLength - 50);

  for (const OptionDescriptionPtr& odPtr : sortedOptions) {
    printedDescriptions.add(odPtr);
  }

  os << this->usagePrefix() << "\n\n" << printedDescriptions << std::endl;
  return os.str();
}

std::string Options::usage() const
{
  OptionDescriptionPtrs sortedOptions = namedOptions_.options();
  std::sort(sortedOptions.begin(), sortedOptions.end(), compareOptionName);

  if (vm_.count(HELP_DEFAULTS_STR)) {
    return helpDefaults(sortedOptions);
  }

  return help(sortedOptions);
}

}  // namespace common
}  // namespace dragenos
