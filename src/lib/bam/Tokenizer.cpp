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

#include "bam/Tokenizer.hpp"
#include "common/Debug.hpp"

namespace dragenos {
namespace bam {

bool Tokenizer::next()
{
  while (true) {
    currentToken_.envelop(&*bufferIterator_);

    const std::size_t bufferLeft = std::distance(bufferIterator_, buffer_.end());
    if (sizeof(BamRecordHeader) > bufferLeft || currentToken_.size() > bufferLeft)  // incomplete
    {
      if (!input_ && !input_.eof()) {
        throw std::ios_base::failure(strerror(errno));
      }

      if (!input_.eof()) {
        const std::size_t available =
            std::distance(buffer_.begin(), bufferIterator_) + buffer_.capacity() - buffer_.size();
        if (!available) {
          throw std::logic_error(
              std::string("Insufficient buffer capacity ")
              << buffer_.capacity() << " to load complete record from stream around offset "
              << input_.tellg());
        }

        buffer_.erase(buffer_.begin(), bufferIterator_);
        bufferIterator_ = buffer_.end();
        buffer_.resize(buffer_.capacity());
        input_.read(&*bufferIterator_, available);
        buffer_.resize(buffer_.capacity() - available + input_.gcount());
        bufferIterator_ = buffer_.begin();

        // reset token before having a chance to throw an exception to avoid invalid iterators
        currentToken_.envelop(&*bufferIterator_);
        const bool complete =
            sizeof(BamRecordHeader) <= buffer_.size() && currentToken_.size() <= buffer_.size();
        if (!input_ && !input_.eof()) {
          throw std::ios_base::failure(strerror(errno));
        }

        // now that we know IO was successful, do the checks on the token
        if (!complete) {
          if (input_.eof()) {
            if (buffer_.end() != bufferIterator_) {
              throw std::logic_error(
                  std::string("Invalid bam record at the end of the stream around offset ")
                  << input_.tellg());
            }
            return false;
          } else {
            throw std::logic_error(
                std::string("Failed to read complete record into buffer of capacity ")
                << buffer_.capacity() << " around offset " << input_.tellg() << " Buffer too small?");
          }
        }
      } else if (buffer_.end() != bufferIterator_) {
        throw std::logic_error(
            std::string("Invalid bam record at the end of the stream around offset ") << input_.tellg());
      } else {
        return false;
      }
    }

    //  std::cerr << "complete:" << currentToken_ << std::endl;
    bufferIterator_ += currentToken_.size();
    if (currentToken_.secondary() || currentToken_.suplementary()) {
      // skip the rubbish
      continue;
    }
    return true;
  }
}

}  // namespace bam
}  // namespace dragenos
