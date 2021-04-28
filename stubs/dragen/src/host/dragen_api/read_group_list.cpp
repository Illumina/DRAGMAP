#include "read_group_list.hpp"

ReadGroupList& GetReadGroupList()
{
  static ReadGroupList readGroupList;
  return readGroupList;
}
