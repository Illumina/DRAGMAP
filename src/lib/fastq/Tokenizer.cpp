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

#include "fastq/Tokenizer.hpp"
#include "common/Debug.hpp"

namespace dragenos {
namespace fastq {

bool Tokenizer::next()
{
  if (!currentToken_.reset(bufferIterator_, buffer_.end()))  // incomplete
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
            << buffer_.capacity() << " to load complete record from stream around offset " << input_.tellg()
            << " token:" << currentToken_);
      }
      //      std::cerr << "erase for " << std::distance(bufferIterator_, buffer_.end()) << "\n";
      //      bufferIterator_ = buffer_.erase(buffer_.begin(), bufferIterator_);
      //      buffer_.resize(buffer_.capacity());
      //      input_.read(&*(buffer_.end() - available), available);
      //      buffer_.resize(buffer_.capacity() - available + input_.gcount());

      // this potentially spends less time zeroing-out bytes though no evidence seen
      std::move(bufferIterator_, buffer_.end(), buffer_.begin());
      bufferIterator_ = buffer_.begin() + std::distance(bufferIterator_, buffer_.end());
      buffer_.resize(buffer_.capacity());
      input_.read(&*bufferIterator_, available);
      buffer_.resize(buffer_.capacity() - available + input_.gcount());
      bufferIterator_ = buffer_.begin();

      if (mixedNewline_) {
        std::replace(buffer_.begin() + buffer_.capacity() - available, buffer_.end(), '\r', '\n');
      }

      // reset token before having a chance to throw an exception to avoid invalid iterators
      const bool complete = currentToken_.reset(bufferIterator_, buffer_.end());
      if (!input_ && !input_.eof()) {
        throw std::ios_base::failure(strerror(errno));
      }

      // now that we know IO was successful, do the checks on the token
      if (complete) {
        // move on
        bufferIterator_ = currentToken_.end();
      } else if (input_.eof()) {
        if (!currentToken_.empty()) {
          throw std::logic_error(
              std::string("Invalid fastq record at the end of the stream around offset ")
              << input_.tellg() << " token:" << currentToken_);
        }
      } else {
        throw std::logic_error(
            std::string("Failed to read complete record into buffer of capacity ")
            << buffer_.capacity() << " around offset " << input_.tellg() << " token:" << currentToken_
            << " Buffer too small?");
      }
    } else if (!currentToken_.empty()) {
      throw std::logic_error(
          std::string("Invalid fastq record at the end of the stream around offset ")
          << input_.tellg() << " token:" << currentToken_);
    } else {
      assert(!currentToken_.valid());
    }
  } else {
    bufferIterator_ = currentToken_.end();
  }
  return currentToken_.valid();
}

}  // namespace fastq
}  // namespace dragenos
