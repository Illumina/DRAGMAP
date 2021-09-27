#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/program_options.hpp>

class Sam {
public:
  Sam &operator=(const std::string &line);
  bool operator==(const Sam &rhs) const;
  const std::string &getUnparsed() const { return unparsed_; }
  const std::string &getName() const { return name_; }
  unsigned getFlags() const { return flags_; }
  unsigned getFlag(unsigned i) const { return flags_ & (1 << i); }
  const std::string &getSequence() const { return sequence_; }
  long long getPosition() const { return position_; }
  int getMapq() const { return mapq_; }
  const std::string &getCigar() const { return cigar_; }
  const std::unordered_map<std::string, std::string> &getTags() const {
    return tags_;
  }
  // reset
  void reset() { tags_.clear(); }
  // all individual flags
  bool hasMultipleSegments() const { return flags_ & 0x1; }
  bool allSegmentsAreAligned() const { return flags_ & 0x2; }
  bool isUnmapped() const { return flags_ & 0x4; }
  bool nextIsUnmapped() const { return flags_ & 0x8; }
  bool isReverveComplemented() const { return flags_ & 0x10; }
  bool nextIsReverveComplemented() const { return flags_ & 0x20; }
  bool isFirst() const { return flags_ & 0x40; }
  bool isLast() const { return flags_ & 0x80; }
  bool isSecondary() const { return flags_ & 0x100; }
  bool isFiltered() const { return flags_ & 0x200; }
  bool isDuplicate() const { return flags_ & 0x400; }
  bool isSupplementary() const { return flags_ & 0x800; }

private:
  std::string unparsed_;
  std::string name_;
  unsigned flags_;
  std::string sequence_;
  long long position_;
  std::string cigar_;
  std::unordered_map<std::string, std::string> tags_;
  int mapq_;
};

typedef std::pair<unsigned, unsigned> Count;
typedef std::pair<std::string, unsigned>
    AltID; // for secondary/supplementary alignments

// tmp containers for secondary/supplementary alignments
typedef std::unordered_map<AltID, Sam, boost::hash<AltID>> SamBuffer;

struct Statistics {
  Statistics() {
    position.fill(Count(0, 0));
    flags.fill(Count(0, 0));
    mapq.fill(Count(0, 0));
    mapqDirection.fill(Count(0, 0));
    cigar.fill(Count(0, 0));
  }
  void write(const std::string outputDirectory) const;
  static unsigned mapqBin(const int mapq);
  unsigned missing = 0;
  unsigned extra = 0;
  unsigned match = 0;
  unsigned mismatch = 0;
  unsigned moreAlignments = 0;
  unsigned lessAlignments = 0;
  Count unmapped = {0, 0};
  Count mapped = {0, 0};
  std::array<Count, 4> position; // MAPQ: 60, 30-59, 1-29, 0
  std::array<Count, 12> flags;
  std::array<Count, 4> mapq;
  std::array<Count, 4> mapqDirection;
  std::array<Count, 4> cigar;

  // note that a record may contribute to multiple mismatched tags
  std::unordered_map<std::string, std::array<Count, 4>> tags;
  std::unordered_map<std::string, int> missingTags;
  std::unordered_map<std::string, int> extraTags;

  // secondary alignments stats
  unsigned secondaryMissing = 0;
  unsigned secondaryExtra = 0;
  unsigned secondaryMatch = 0;

  // supplementary alignments stats
  unsigned supplementaryMissing = 0;
  unsigned supplementaryExtra = 0;
  unsigned supplementaryMatch = 0;
};

bool getSamRecords(std::istream &is, std::vector<Sam> &samRecords,
                   SamBuffer &secondaryBuffers, SamBuffer &supplementaryBuffers,
                   Sam &next);
void compare(
    std::vector<Sam> &referenceRecords, std::vector<Sam> &newRecords,
    Statistics &statistics, std::ostream &missingMappings,
    std::ostream &extraMappings, std::ostream &positionMismatch,
    std::ostream &positionMismatchReference, std::ostream &flagMismatch,
    std::ostream &flagMismatchReference, std::ostream &cigarMismatch,
    std::ostream &cigarMismatchReference, std::ostream &pairEndInfoMismatch,
    std::ostream &pairEndInfoMismatchReference,
    std::map<std::string, std::shared_ptr<std::ostream>> &tagMismatch,
    std::map<std::string, std::shared_ptr<std::ostream>> &tagMismatchReference);
std::ostream &operator<<(std::ostream &os, const Sam &sam);
std::ostream &operator<<(std::ostream &os, const Statistics &statistics);

void compareSecondary(SamBuffer &refSamBuffer, SamBuffer &newSamBuffer,
                      Statistics &statistics);
void compareSupplementary(SamBuffer &refSamBuffer, SamBuffer &newSamBuffer,
                          Statistics &statistics);

void makeMismatchFiles(
    const std::string &tag, const std::string &fileSuffix,
    const boost::filesystem::path &outputDirectory,
    std::map<std::string, std::shared_ptr<std::ostream>> &tagMismatch) {
  tagMismatch.emplace(
      tag + std::to_string(Statistics::mapqBin(0)),
      std::make_shared<std::ofstream>(
          (outputDirectory / ("tag" + tag + "Mismatch0" + fileSuffix))
              .string()));
  tagMismatch.emplace(
      tag + std::to_string(Statistics::mapqBin(1)),
      std::make_shared<std::ofstream>(
          (outputDirectory / ("tag" + tag + "Mismatch1-29" + fileSuffix))
              .string()));
  tagMismatch.emplace(
      tag + std::to_string(Statistics::mapqBin(30)),
      std::make_shared<std::ofstream>(
          (outputDirectory / ("tag" + tag + "Mismatch30-59" + fileSuffix))
              .string()));
  tagMismatch.emplace(
      tag + std::to_string(Statistics::mapqBin(60)),
      std::make_shared<std::ofstream>(
          (outputDirectory / ("tag" + tag + "Mismatch60" + fileSuffix))
              .string()));
}

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;
  po::options_description options(
      "Usage: compare reference.sam new.sam [options]\n"
      "compare the specified 'new.sam' file against 'reference.sam'.\n"
      "Supported options are:");
  options.add_options()("help", "produce help message")(
      "key-value,k", po::value<bool>()->default_value(false),
      "specify the output should be a list of key=value on stdout instead of a "
      "set of cvs files")("output,o",
                          po::value<std::string>()->default_value("./"),
                          "specify the output directory")(
      "position-mismatches", po::value<bool>()->default_value(false),
      "specifies that the list of records with MAPQ>0 and incorrect positions "
      "should be produced")(
      "cigar-mismatches", po::value<bool>()->default_value(false),
      "specifies that the list of records with MAPQ>0 and correct positions "
      "and mismatching CIGAR strings should be produced")(
      "flag-mismatches", po::value<bool>()->default_value(false),
      "specifies that the list of records with inconsistent-flag mappings "
      "should be produced")(
      "missing-mappings", po::value<bool>()->default_value(false),
      "specifies that the list of records with missing mappings should be "
      "produced")("extra-mappings", po::value<bool>()->default_value(false),
                  "specifies that the list of records with extra mappings "
                  "should be produced")(
      "pair-end-mismatches", po::value<bool>()->default_value(false),
      "speicifes that the list of records that are inconsistent from pair-end "
      "mapping(sam files are expected to be sorted by read name and FR order)")(
      "input-file,i", po::value<std::vector<std::string>>(),
      "input SAM files reference.sam and new.sam -- these can also be "
      "specified simply as positional options");
  po::positional_options_description positional;
  positional.add("input-file", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
                .options(options)
                .positional(positional)
                .run(),
            vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << options << "\n";
    return 1;
  }
  const std::vector<std::string> inputFiles =
      vm["input-file"].as<std::vector<std::string>>();
  if (2 != inputFiles.size()) {
    std::cerr << "ERROR: exactly two input files must be specified: "
                 "reference.sam and new.sam"
              << std::endl;
    std::cout << options << "\n";
    return 1;
  }
  std::array<std::ifstream, 2> is;
  for (auto i : {0, 1}) {
    is[i].open(argv[i + 1]);
    if (!is[i]) {
      std::cerr << "failed to open " << inputFiles[i] << ": " << errno << ": "
                << strerror(errno) << std::endl;
      exit(2);
    }
  }

  std::array<SamBuffer, 2> secondaryBuffers;
  std::array<SamBuffer, 2> supplementaryBuffers;

  // skip the headers
  std::array<Sam, 2> nextRecords;
  for (auto i : {0, 1}) {
    std::string line;
    while (getline(is[i], line) && ('@' == line[0])) {
    }
    nextRecords[i] = line;
  }
  // compare
  std::array<std::vector<Sam>, 2> samRecords;

  getSamRecords(is[0], samRecords[0], secondaryBuffers[0],
                supplementaryBuffers[0], nextRecords[0]);
  getSamRecords(is[1], samRecords[1], secondaryBuffers[1],
                supplementaryBuffers[1], nextRecords[1]);
  Statistics statistics;
  // while (getSamRecords(is[0], samRecords[0], nextRecords[0]) ||
  // getSamRecords(is[1], samRecords[1], nextRecords[1]))
  std::ofstream positionMismatch;
  std::ofstream positionMismatchReference;
  const boost::filesystem::path outputDirectory =
      vm["output"].as<std::string>();
  if (vm["position-mismatches"].as<bool>()) {
    positionMismatch.open((outputDirectory / "positionMismatch.sam").c_str());
    positionMismatchReference.open(
        (outputDirectory / "positionMismatchReference.sam").c_str());
  }
  std::ofstream cigarMismatch;
  std::ofstream cigarMismatchReference;
  if (vm["cigar-mismatches"].as<bool>()) {
    cigarMismatch.open((outputDirectory / "cigarMismatch.sam").c_str());
    cigarMismatchReference.open(
        (outputDirectory / "cigarMismatchReference.sam").c_str());
  }
  std::ofstream missingMappings;
  if (vm["missing-mappings"].as<bool>()) {
    missingMappings.open((outputDirectory / "missingMappings.sam").c_str());
  }

  std::ofstream extraMappings;
  if (vm["extra-mappings"].as<bool>()) {
    extraMappings.open((outputDirectory / "extraMappings.sam").c_str());
  }

  std::ofstream flagMismatch;
  std::ofstream flagMismatchReference;
  if (vm["flag-mismatches"].as<bool>()) {
    flagMismatch.open((outputDirectory / "flagMismatch.sam").c_str());
    flagMismatchReference.open(
        (outputDirectory / "flagMismatchReference.sam").c_str());
  }

  std::ofstream pairEndInfoMismatch;
  std::ofstream pairEndInfoMismatchReference;
  if (vm["pair-end-mismatches"].as<bool>()) {
    pairEndInfoMismatch.open((outputDirectory / "pairEndMismatch.sam").c_str());
    pairEndInfoMismatchReference.open(
        (outputDirectory / "pairEndMismatchReference.sam").c_str());
  }

  std::map<std::string, std::shared_ptr<std::ostream>> tagMismatch;
  std::map<std::string, std::shared_ptr<std::ostream>> tagMismatchReference;
  makeMismatchFiles("XS", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("XS", "Reference.sam", outputDirectory,
                    tagMismatchReference);
  makeMismatchFiles("NM", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("NM", "Reference.sam", outputDirectory,
                    tagMismatchReference);
  makeMismatchFiles("AS", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("AS", "Reference.sam", outputDirectory,
                    tagMismatchReference);
  makeMismatchFiles("XQ", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("XQ", "Reference.sam", outputDirectory,
                    tagMismatchReference);
  makeMismatchFiles("SA", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("SA", "Reference.sam", outputDirectory,
                    tagMismatchReference);
  makeMismatchFiles("MAPQ", ".sam", outputDirectory, tagMismatch);
  makeMismatchFiles("MAPQ", "Reference.sam", outputDirectory,
                    tagMismatchReference);

  std::vector<Sam> empty;
  int counter = 0;
  while ((!samRecords[0].empty()) || (!samRecords[1].empty())) {
    // make sure that the comparison is for the same read name of against an
    // empty set
    if (counter % 1000000 == 0)
      std::cerr << "Input " << counter << " Records(Single or Pair)...\n";
    counter++;
    compareSecondary(secondaryBuffers[0], secondaryBuffers[1], statistics);
    compareSupplementary(supplementaryBuffers[0], supplementaryBuffers[1],
                         statistics);
    if (samRecords[0].empty() || samRecords[1].empty() ||
        (samRecords[0][0].getName() == samRecords[1][0].getName())) {
      compare(samRecords[0], samRecords[1], statistics, missingMappings,
              extraMappings, positionMismatch, positionMismatchReference,
              flagMismatch, flagMismatchReference, cigarMismatch,
              cigarMismatchReference, pairEndInfoMismatch,
              pairEndInfoMismatchReference, tagMismatch, tagMismatchReference);
      getSamRecords(is[0], samRecords[0], secondaryBuffers[0],
                    supplementaryBuffers[0], nextRecords[0]);
      getSamRecords(is[1], samRecords[1], secondaryBuffers[1],
                    supplementaryBuffers[1], nextRecords[1]);
    } else if (samRecords[0][0].getName() < samRecords[1][0].getName()) {
      compare(samRecords[0], empty, statistics, missingMappings, extraMappings,
              positionMismatch, positionMismatchReference, flagMismatch,
              flagMismatchReference, cigarMismatch, cigarMismatchReference,
              pairEndInfoMismatch, pairEndInfoMismatchReference, tagMismatch,
              tagMismatchReference);
      getSamRecords(is[0], samRecords[0], secondaryBuffers[0],
                    supplementaryBuffers[0], nextRecords[0]);
    } else // (samRecords[0].getName() > samRecords[1].getName())
    {
      compare(empty, samRecords[1], statistics, missingMappings, extraMappings,
              positionMismatch, positionMismatchReference, flagMismatch,
              flagMismatchReference, cigarMismatch, cigarMismatchReference,
              pairEndInfoMismatch, pairEndInfoMismatchReference, tagMismatch,
              tagMismatchReference);
      getSamRecords(is[1], samRecords[1], secondaryBuffers[1],
                    supplementaryBuffers[1], nextRecords[1]);
    }
  }
  if (vm["key-value"].as<bool>()) {
    std::cout << statistics << std::endl;
  } else {
    statistics.write(vm["output"].as<std::string>());
  }
}

Sam &Sam::operator=(const std::string &line) {
  reset();
  unparsed_ = line;
  /*
  std::istringstream is(line);
  is >> name_;
  is >> flags_;
  is >> sequence_;
  is >> position_;
  is >> mapq_;
  is >> cigar_;
  std::string dummy;
  is >> dummy; // RNEXT
  is >> dummy; // PNEXT
  is >> dummy; // TLEN
  is >> dummy; // SEQ
  is >> dummy; // QUAL
  */
  std::vector<std::string> items;
  std::vector<std::string> results;
  boost::algorithm::split(items, line, boost::is_any_of("\t"));
  auto iter = items.begin();
  name_ = *iter++;
  flags_ = std::stoi(*iter++);
  sequence_ = *iter++;
  position_ = std::stoi(*iter++);
  mapq_ = std::stoi(*iter++);
  cigar_ = *iter++;
  iter += 5;
  std::string dummy;
  while (iter != items.end()) {
    dummy = *iter++;
    boost::algorithm::split(results, dummy, boost::is_any_of(":"));
    tags_[results[0]] = results[2];
  }
  return *this;
}

bool Sam::operator==(const Sam &rhs) const {
  std::unordered_set<std::string> tagsBlockList = {"XQ", "RG"};
  // if(tags_.size() != rhs.getTags().size()) return false;
  for (const auto &kv : tags_) {
    if (tagsBlockList.find(kv.first) != tagsBlockList.end())
      continue;
    if (rhs.getTags().find(kv.first) != rhs.getTags().end() and
        tags_.at(kv.first) != rhs.getTags().at(kv.first))
      return false;
  }

  return name_ == rhs.getName() and flags_ == rhs.getFlags() and
         sequence_ == rhs.getSequence() and position_ == rhs.getPosition() and
         mapq_ == rhs.getMapq() and cigar_ == rhs.getCigar();
}

bool getSamRecords(std::istream &is, std::vector<Sam> &samRecords,
                   SamBuffer &secondaryBuffers, SamBuffer &supplementaryBuffers,
                   Sam &next) {
  samRecords.clear();
  secondaryBuffers.clear();
  supplementaryBuffers.clear();
  if (!is) {
    return false;
  }
  if (next.isSecondary()) {
    secondaryBuffers.insert(std::make_pair(
        std::make_pair(next.getName() + (next.isFirst() ? "1" : "2") +
                           next.getCigar() + next.getSequence(),
                       next.getPosition()),
        next));
  } else if (next.isSupplementary()) {
    supplementaryBuffers.insert(std::make_pair(
        std::make_pair(next.getName() + (next.isFirst() ? "1" : "2") +
                           next.getCigar() + next.getSequence(),
                       next.getPosition()),
        next));
  } else {
    samRecords.push_back(next);
  }
  std::string line;
  Sam last = next;
  while (getline(is, line) && ((next = line).getName() == last.getName())) {
    if (next.isSecondary()) {
      secondaryBuffers.insert(std::make_pair(
          std::make_pair(next.getName() + (next.isFirst() ? "1" : "2") +
                             next.getCigar() + next.getSequence(),
                         next.getPosition()),
          next));
    } else if (next.isSupplementary()) {
      supplementaryBuffers.insert(std::make_pair(
          std::make_pair(next.getName() + (next.isFirst() ? "1" : "2") +
                             next.getCigar() + next.getSequence(),
                         next.getPosition()),
          next));
    } else {
      samRecords.push_back(next);
    }
    last = next;
  }
  return true;
}

void compareSingle(
    const Sam &refRecord, const Sam &newRecord, Statistics &statistics,
    std::ostream &positionMismatch, std::ostream &positionMismatchReference,
    std::ostream &flagMismatch, std::ostream &flagMismatchReference,
    std::ostream &cigarMismatch, std::ostream &cigarMismatchReference,
    std::map<std::string, std::shared_ptr<std::ostream>> &tagMismatch,
    std::map<std::string, std::shared_ptr<std::ostream>>
        &tagMismatchReference) {
  if (refRecord.isUnmapped()) {
    if (newRecord.isUnmapped()) {
      ++statistics.unmapped.first;
    } else {
      ++statistics.unmapped.second;
    }
    return;
  } else {
    if (!newRecord.isUnmapped()) {
      ++statistics.mapped.first;
      // both mapped, keep comparing
    } else {
      ++statistics.mapped.second;
      return;
    }
  }

  if (refRecord == newRecord)
    ++statistics.match;
  else
    ++statistics.mismatch;

  const unsigned mapqBin = Statistics::mapqBin(refRecord.getMapq());
  // if ((referenceRecords[0].getSequence() == newRecords[0].getSequence()) &&
  // (referenceRecords[0].getPosition() == newRecords[0].getPosition()))
  if ((refRecord.getSequence() == newRecord.getSequence()) &&
      (refRecord.getPosition() + 200 > newRecord.getPosition()) &&
      (refRecord.getPosition() < newRecord.getPosition() + 200)) {
    ++statistics.position[mapqBin].first;
    // both mapped at the same location - keep comparing
  } else {
    ++statistics.position[mapqBin].second;
    if (positionMismatch && (3 != mapqBin)) {
      positionMismatch << newRecord.getUnparsed() << std::endl;
      positionMismatchReference << refRecord.getUnparsed() << std::endl;
    }
    return;
  }
  for (unsigned i = 0; i < 12; ++i) {
    if (refRecord.getFlag(i) == newRecord.getFlag(i)) {
      ++statistics.flags[i].first;
    } else {
      ++statistics.flags[i].second;
      if (flagMismatch && (3 != mapqBin)) {
        flagMismatch << newRecord.getUnparsed() << std::endl;
        flagMismatchReference << refRecord.getUnparsed() << std::endl;
      }
    }
  }
  if (refRecord.getMapq() == newRecord.getMapq()) {
    ++statistics.mapq[mapqBin].first;
  } else {
    ++statistics.mapq[mapqBin].second;
    if (refRecord.getMapq() < newRecord.getMapq()) {
      ++statistics.mapqDirection[mapqBin].first;
    } else {
      ++statistics.mapqDirection[mapqBin].second;
    }
    const std::string key = "MAPQ" + std::to_string(mapqBin);
    const auto stream = tagMismatch.find(key);
    if (tagMismatch.end() != stream) {
      *stream->second << newRecord.getUnparsed() << std::endl;
      *tagMismatchReference.at(key) << refRecord.getUnparsed() << std::endl;
    }
  }
  if (refRecord.getCigar() == newRecord.getCigar()) {
    ++statistics.cigar[mapqBin].first;
  } else {
    ++statistics.cigar[mapqBin].second;
    if (cigarMismatch && (3 != mapqBin)) {
      cigarMismatch << newRecord.getUnparsed() << std::endl;
      cigarMismatchReference << refRecord.getUnparsed() << std::endl;
    }
  }

  // compare tags
  if (not refRecord.getTags().empty()) {
    for (const auto &kv : refRecord.getTags()) {
      if (newRecord.getTags().find(kv.first) != newRecord.getTags().end()) {
        if (newRecord.getTags().at(kv.first) ==
            refRecord.getTags().at(kv.first)) {
          statistics.tags[kv.first].at(mapqBin).first++;
        } else {
          statistics.tags[kv.first].at(mapqBin).second++;
          const std::string key = kv.first + std::to_string(mapqBin);
          const auto stream = tagMismatch.find(key);
          if (tagMismatch.end() != stream) {
            *stream->second << newRecord.getUnparsed() << std::endl;
            *tagMismatchReference.at(key)
                << refRecord.getUnparsed() << std::endl;
          }
        }
      } else {
        statistics.missingTags[kv.first]++;
      }
    }
    for (const auto &kv : newRecord.getTags()) {
      if (refRecord.getTags().find(kv.first) == refRecord.getTags().end()) {
        statistics.extraTags[kv.first]++;
      }
    }
  } else {
    for (const auto &kv : newRecord.getTags()) {
      statistics.extraTags[kv.first]++;
    }
  }
}

void compare(std::vector<Sam> &referenceRecords, std::vector<Sam> &newRecords,
             Statistics &statistics, std::ostream &missingMappings,
             std::ostream &extraMappings, std::ostream &positionMismatch,
             std::ostream &positionMismatchReference,
             std::ostream &flagMismatch, std::ostream &flagMismatchReference,
             std::ostream &cigarMismatch, std::ostream &cigarMismatchReference,
             std::ostream &pairEndInfoMismatch,
             std::ostream &pairEndInfoMismatchReference,
             std::map<std::string, std::shared_ptr<std::ostream>> &tagMismatch,
             std::map<std::string, std::shared_ptr<std::ostream>>
                 &tagMismatchReference) {
  if (referenceRecords.empty() && newRecords.empty()) {
    std::cerr << "No SAM record to compare" << std::endl;
    exit(3);
  }
  if (newRecords.empty()) {
    for (const auto &refRec : referenceRecords) {
      if (missingMappings) {
        missingMappings << refRec.getUnparsed() << std::endl;
      }
      ++statistics.missing;
    }
    return;
  }
  if (referenceRecords.empty()) {
    for (const auto &newRec : newRecords) {
      if (extraMappings) {
        extraMappings << newRec.getUnparsed() << std::endl;
      }
      ++statistics.extra;
    }
    return;
  }

  if (referenceRecords.size() == 2 &&
      newRecords.size() ==
          2) // TODO: append additional pair-end statistics if needed
  {
    // if(!referenceRecords[0].isFirst()) std::swap(referenceRecords[0],
    // referenceRecords[1]); if(!newRecords[0].isFirst())
    // std::swap(newRecords[0], newRecords[1]);

    // if(referenceRecords[0].getInsertSize() != newRecords[0].getInsertSize())
    // {
    //   pairEndInfoMismatch << newRecords[0].getUnparsed() << std::endl;
    //   pairEndInfoMismatch << newRecords[1].getUnparsed() << std::endl;
    //   pairEndInfoMismatchReference << referenceRecords[0].getUnparsed() <<
    //   std::endl; pairEndInfoMismatchReference <<
    //   referenceRecords[1].getUnparsed() << std::endl;
    // }
  } else if (referenceRecords.size() == 1 && newRecords.size() == 1) {
    if (referenceRecords[0].getFlag(6) != newRecords[0].getFlag(6)) {
      // if (missingMappings)
      // {
      //   missingMappings << referenceRecords[0].getUnparsed() << std::endl;
      // }
      // ++statistics.missing;
      // ++statistics.extra;
      std::cerr << "ERROR: Flags of pair-end reads detected, but "
                   "referenceRecords file only has one end."
                << std::endl;
      exit(3);
    }
  } else if (referenceRecords.size() == 2 && newRecords.size() == 1) {
    Sam missingRec = referenceRecords[0].getFlag(6) == newRecords[0].getFlag(6)
                         ? referenceRecords[1]
                         : referenceRecords[0];
    if (missingMappings) {
      missingMappings << missingRec.getUnparsed() << std::endl;
    }
    ++statistics.missing;
  } else if (referenceRecords.size() == 1 && newRecords.size() == 2) {
    Sam newRec = referenceRecords[0].getFlag(6) == newRecords[0].getFlag(6)
                     ? newRecords[1]
                     : newRecords[0];
    if (extraMappings) {
      extraMappings << newRec.getUnparsed() << std::endl;
    }
    ++statistics.extra;
  }

  // handle single comparison
  for (const auto &refRec : referenceRecords) {
    for (const auto &newRec : newRecords) {
      if (refRec.getFlag(6) == newRec.getFlag(6)) // both or neither
      {
        compareSingle(refRec, newRec, statistics, positionMismatch,
                      positionMismatchReference, flagMismatch,
                      flagMismatchReference, cigarMismatch,
                      cigarMismatchReference, tagMismatch,
                      tagMismatchReference);
      }
    }
  }
}

void compareSecondary(SamBuffer &refSamBuffer, SamBuffer &newSamBuffer,
                      Statistics &statistics) {
  for (const auto &kv : refSamBuffer) {
    if (newSamBuffer.find(kv.first) != newSamBuffer.end())
      statistics.secondaryMatch++;
    else
      statistics.secondaryMissing++;
  }

  for (const auto &kv : newSamBuffer) {
    if (refSamBuffer.find(kv.first) == refSamBuffer.end())
      statistics.secondaryExtra++;
  }
}

void compareSupplementary(SamBuffer &refSamBuffer, SamBuffer &newSamBuffer,
                          Statistics &statistics) {
  for (const auto &kv : refSamBuffer) {
    if (newSamBuffer.find(kv.first) != newSamBuffer.end())
      statistics.supplementaryMatch++;
    else
      statistics.supplementaryMissing++;
  }

  for (const auto &kv : newSamBuffer) {
    if (refSamBuffer.find(kv.first) == refSamBuffer.end())
      statistics.supplementaryExtra++;
  }
}

std::ostream &operator<<(std::ostream &os, const Statistics &statistics) {
  os << "missing=" << statistics.missing << std::endl;
  os << "extra=" << statistics.extra << std::endl;
  os << "match=" << statistics.match << std::endl;
  os << "mismatch=" << statistics.mismatch << std::endl;
  static const std::array<std::string, 4> mapqStrings{"60", "30-59", "1-29",
                                                      "0"};
  for (auto i : {0, 1, 2, 3}) {
    os << "position_match_" << mapqStrings[i] << "="
       << statistics.position[i].first << std::endl;
    os << "position_mismatch_" << mapqStrings[i] << "="
       << statistics.position[i].second << std::endl;
  }
  for (unsigned i = 0; 12 > i; ++i) {
    os << "flag" << i << "_match_0-60"
       << "=" << statistics.flags[i].first << std::endl;
    os << "flag" << i << "_mismatch_0-60"
       << "=" << statistics.flags[i].second << std::endl;
  }
  for (auto i : {0, 1, 2, 3}) {
    os << "mapq_match_" << mapqStrings[i] << "=" << statistics.mapq[i].first
       << std::endl;
    os << "mapq_mismatch_" << mapqStrings[i] << "=" << statistics.mapq[i].second
       << std::endl;
    os << "mapq_HW_LT_SW_" << mapqStrings[i] << "="
       << statistics.mapqDirection[i].first << std::endl;
    os << "mapq_HW_GT_SW_" << mapqStrings[i] << "="
       << statistics.mapqDirection[i].second << std::endl;
    os << "mapqx_match_" << mapqStrings[i] << "= NOT IMPLEMENTED" << std::endl;
    os << "mapqx_mismatch_" << mapqStrings[i] << "= NOT IMPLEMENTED"
       << std::endl;
  }
  for (auto i : {0, 1, 2, 3}) {
    os << "cigar_match_" << mapqStrings[i] << "=" << statistics.cigar[i].first
       << std::endl;
    os << "cigar_mismatch_" << mapqStrings[i] << "="
       << statistics.cigar[i].second << std::endl;
  }

  for (auto i : {0, 1, 2, 3}) {
    for (const auto &kv : statistics.tags) {
      os << "tag_" << kv.first << "_match_" << mapqStrings[i] << "="
         << kv.second.at(i).first << std::endl;
      os << "tag_" << kv.first << "_mismatch_" << mapqStrings[i] << "="
         << kv.second.at(i).second << std::endl;
    }
  }

  for (const auto &kv : statistics.extraTags) {
    os << "tag_" << kv.first << "_extra"
       << "=" << kv.second << std::endl;
  }

  for (const auto &kv : statistics.missingTags) {
    os << "tag_" << kv.first << "_missing_"
       << "=" << kv.second << std::endl;
  }

  os << "supplementary_match=" << statistics.supplementaryMatch << std::endl;
  os << "supplementary_missing=" << statistics.supplementaryMissing
     << std::endl;
  os << "supplementary_extra=" << statistics.supplementaryExtra << std::endl;
  os << "secondary_match=" << statistics.secondaryMatch << std::endl;
  os << "secondary_missing=" << statistics.secondaryMissing << std::endl;
  os << "secondary_extra=" << statistics.secondaryExtra << std::endl;
  return os;
}

unsigned Statistics::mapqBin(const int mapq) {
  if (255 == mapq)
    return 3; // MAPQ not available
  if (mapq >= 60)
    return 0;
  if (mapq >= 30)
    return 1;
  if (mapq > 0)
    return 2;
  return 3;
}

std::ostream &operator<<(std::ostream &os, const Sam &sam) {
  return os << sam.getName() << '\t' << sam.getFlags() << '\t' << std::endl;
}

void Statistics::write(const std::string outputDirectory) const {
  std::cerr << "writing into " << outputDirectory << std::endl;
  const boost::filesystem::path path = outputDirectory;
  {
    std::ofstream os((path / "unmapped.csv").c_str());
    if (!(os << "missing,extra\n" << missing << "," << extra << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  {
    std::ofstream os((path / "total_match_mismatch.csv").c_str());
    if (!(os << "match,mismatch\n" << match << "," << mismatch << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  const std::vector<std::string> mapqBins{"60", "30_59", "1_29", "0"};
  for (unsigned i = 0; 4 > i; ++i) {
    const std::string label = std::string("position_mismatch_") + mapqBins[i];
    std::ofstream os((path / (label + ".csv")).c_str());
    if (!(os << label << "\n" << position[i].second << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  for (unsigned i = 0; 4 > i; ++i) {
    const std::string label = std::string("position_match_") + mapqBins[i];
    std::ofstream os((path / (label + ".csv")).c_str());
    if (!(os << label << "\n" << position[i].first << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  {
    std::ofstream os((path / "flag_mismatch.csv").c_str());
    os << "flag_" << 0;
    for (unsigned i = 1; 12 > i; ++i) {
      os << ",flag_" << i;
    }
    os << std::endl << flags[0].second;
    for (unsigned i = 1; 12 > i; ++i) {
      os << "," << flags[i].second;
    }
    if (!(os << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  for (unsigned i = 0; 4 > i; ++i) {
    const std::string label = std::string("mapq_mismatch_") + mapqBins[i];
    std::ofstream os((path / (label + ".csv")).c_str());
    if (!(os << label << ",HW_LT_SW_" << mapqBins[i] << ",HW_GT_SW_"
             << mapqBins[i] << "\n"
             << mapq[i].second << "," << mapqDirection[i].first << ","
             << mapqDirection[i].second << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
  for (unsigned i = 0; 4 > i; ++i) {
    const std::string label = std::string("cigar_mismatch_") + mapqBins[i];
    std::ofstream os((path / (label + ".csv")).c_str());
    if (!(os << label << "\n" << cigar[i].second << std::endl)) {
      std::cerr << "failed to write output csv file: " << strerror(errno)
                << std::endl;
      exit(5);
    }
  }
}
