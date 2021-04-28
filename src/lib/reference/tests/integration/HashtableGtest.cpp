#include "gtest/gtest.h"

#include "reference/Hashtable.hpp"

/**
 ** \brief integration tests for the hashtable
 **
 ** Note: most of the meaningful tests are currently into the systems tests.
 ** TODO: figure out a way to configure and initialize a hashtable without
 ** relying on external files
 **/
TEST(Hashtable, Constants)
{
  using dragenos::reference::HashRecord;
  using dragenos::reference::Hashtable;
  ASSERT_EQ(8, Hashtable::getRecordsPerBucket());
  ASSERT_EQ(8, Hashtable::getBytesPerRecord());
  ASSERT_EQ(sizeof(HashRecord), Hashtable::getBytesPerRecord());
  ASSERT_EQ(64, Hashtable::getBytesPerBucket());
}
