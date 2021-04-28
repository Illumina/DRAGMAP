// Copyright 2015 Edico Genome Corporation. All rights reserved.
//
// This file contains confidential and proprietary information of the Edico
// Genome Corporation and is protected under the U.S. and international
// copyright and other intellectual property laws.
//
// $Id$
// $Author$
// $Change$
// $DateTime$
//

#ifndef __READ_GROUP_LIST_HPP__
#define __READ_GROUP_LIST_HPP__

#include "dragen_exception.hpp"
class ReadGroupList {
public:
  const std::string &getReadGroupName(const uint16_t idx) const {
    static std::string all = "all";
    return all;
  }
};

// Poor-man's Singleton
ReadGroupList &GetReadGroupList();

#endif
