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

#pragma once

std::string getCurrentExecutableDir();

// Translate config #filePath# to be relative the current user's build system directory.
// Used on jenkins when running the mock DMA suite, where there is only a build directory
// and no /opt/edico at the time the tests execute.
void getBuildConfigPath(std::string& filePath);

void getBuildPath(std::string& filePath, const std::string& subDir);
