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

#pragma once

#include "fastq/FastqNRecordReader.hpp"
#include "options/DragenOsOptions.hpp"

namespace dragenos {
namespace workflow {
void input2Sam(const dragenos::options::DragenOsOptions& options);

}  // namespace workflow
}  // namespace dragenos
