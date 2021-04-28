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

#include "align/Alignment.hpp"
#include "align/InsertSizeParameters.hpp"

namespace dragenos {
namespace align {
int getTlen(
    const dragenos::align::Alignment& r1Align,
    const dragenos::align::Alignment& r2Align,
    InsertSizeParameters::Orientation prevailingOrientation);

}  // namespace align
}  // namespace dragenos
