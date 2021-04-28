// Copyright 2013-2018 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico Genome
// Corporation and is protected under the U.S. and international copyright and other
// intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#include <unistd.h>
#include <limits>
#include <boost/filesystem.hpp>

#include "dragen_exception.hpp"
#include "infra_linux_utils.hpp"
#include "linux_utils.hpp"

//-------------------------------------------------------------------------------swhitmore
// Return the path to the current executable
//
std::string getCurrentExecutableDir()
{
  return infra::getExecutableDirectory();
}

//-------------------------------------------------------------------------------swhitmore
// Translate config #filePath# to be relative the current user's build system directory.
// Used on jenkins when running the mock DMA suite, where there is only a build directory
// and no /opt/edico at the time the tests execute.
//
void getBuildConfigPath(std::string& filePath)
{
  const std::string configDir = "config/";
  getBuildPath(filePath, configDir);
  ASSERT(boost::filesystem::exists(filePath), filePath, " does not exist");
}

// Get the filePath to be relative to the current user's build system directory.
//
void getBuildPath(std::string& filePath, const std::string& subDir)
{
  boost::filesystem::path dragen_exe_dir = infra::getExecutableDirectory();
  boost::filesystem::path p(filePath);
  filePath = dragen_exe_dir.parent_path().string() + "/" + subDir + p.filename().string();
}
