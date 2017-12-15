/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineLz77ViaRlbwt.hpp
 * @brief Online LZ77 computation via online RLBWT.
 * @author Tomohiro I
 * @date 2017-10-12
 */
#ifndef INCLUDE_GUARD_OnlineLz77ViaRlbwt
#define INCLUDE_GUARD_OnlineLz77ViaRlbwt

#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>

// include from Basics
#include "BitsUtil.hpp"
#include "WBitsVec.hpp"
#include "MemUtil.hpp"

// include from BTree
#include "BTree.hpp"
#include "TagRelabelAlgo.hpp"


namespace itmmti
{
  /*!
   * @brief Dynamic run-length encoding supporting access, rank, select, and insert with associated value.
   * @tparam B Parameter of B+tree, which should be in {4, 8, 16, 32, 64, 128}.
   * @par Notation
   *   - T: Current string represented by RLE.
   *   - Mixed tree: B+tree representing RLE of T.
   *     - btmM: Index of bottom node of B+tree on the array-based implementation.
   *             Each bottom node 'btmM' can have 'B' children, which correspond to indexes [btmM * B, (btmM+1) * B).
   *     - idxM: Indexes that are corresponding to children of btmM.
   *   - Separated tree: B+tree separately representing runs for each character.
   *     - btmS: Index of bottom node of B+tree on the array-based implementation (all separated trees share arrays).
   *             Each bottom node 'btmS' can have 'B' children, which correspond to indexes [btmS * B, (btmS+1) * B).
   *     - idxS: Indexes that are corresponding to children of btmS.
   */
  template <uint8_t B = 64, typename AssocT = uint64_t>
  class DynRleAssoc
  {
  public:
    //// Public constant, alias etc.
    using BTreeNodeT = BTreeNode<B>;


  private:
    AssocT * assoc_; //!< Associated value for each leaf.
    typename BTreeNodeT::SuperRootT srootM_; //!< Super root of mixed tree.
    // Information for leaves and elements for mixed tree.
    WBitsVec idxM2S_; //!< Packed array mapping idxM to corresponding idxS.
    BTreeNodeT ** parentM_; //!< Pointer to parent of btmM.
    uint64_t * labelM_; //!< TRA label of btmM.
    uint8_t * idxInSiblingM_; //!< idxInSibling of btmM.
    WBitsVec ** weightVecs_; //!< 'weightVecs_[btmM]' is packed array storing weights of runs under btmM.
    // Alphabet tree: the bottoms of alphabet tree are roots of separated trees.
    typename BTreeNodeT::SuperRootT srootA_; //!< Super root of alphabet tree.
    // Information for leaves and elements for separated tree.
    WBitsVec idxS2M_; //!< Packed array mapping idxS to corresponding idxM.
    BTreeNodeT ** parentS_; //!< Pointer to parent of btmS.
    uint64_t * charS_; //!< 64bit-char of btmS.
    uint8_t * idxInSiblingS_; //!< idxInSibling of btmS.
    uint8_t * numChildrenS_; //!< Num of children of btmS.

    uint8_t traCode_; //!< traCode in [9..16).


  public:
    DynRleAssoc() :
      assoc_(nullptr),
      srootM_(),
      idxM2S_(),
      parentM_(nullptr),
      labelM_(nullptr),
      idxInSiblingM_(nullptr),
      weightVecs_(nullptr),
      srootA_(),
      idxS2M_(),
      parentS_(nullptr),
      charS_(nullptr),
      idxInSiblingS_(nullptr),
      numChildrenS_(nullptr),
      traCode_(9)
    {}


    DynRleAssoc
    (
     const size_t initNumBtms
     ) :
      assoc_(nullptr),
      srootM_(),
      idxM2S_(),
      parentM_(nullptr),
      labelM_(nullptr),
      idxInSiblingM_(nullptr),
      weightVecs_(nullptr),
      srootA_(),
      idxS2M_(),
      parentS_(nullptr),
      charS_(nullptr),
      idxInSiblingS_(nullptr),
      numChildrenS_(nullptr),
      traCode_(9)
    {
      init(initNumBtms);
    }


    ~DynRleAssoc() {
      clearAll();
    }


    /*!
     * @brief Reserve space to accomodate 'initNumBtms' bottoms, and init.
     */
    void init(const size_t initNumBtms) {
      assert(initNumBtms > 0);

      if (isReady()) {
        clearAll();
      }
      reserveBtms(initNumBtms);

      srootM_.setRoot(new BTreeNodeT(reinterpret_cast<BTreeNodeT *>(0), true, true, true, true));
      // sentinel
      parentM_[0] = srootM_.root_;
      idxInSiblingM_[0] = 0;
      labelM_[0] = 0;
      idxM2S_.resize(B);
      idxM2S_.write(0, 0); // sentinel: should not be used
      weightVecs_[0] = new WBitsVec(8, B);
      weightVecs_[0]->resize(1);
      weightVecs_[0]->write(0, 0);
      srootM_.root_->putFirstBtm(reinterpret_cast<BTreeNodeT *>(0), 0);

      // isRoot = true, isBorder = true, isJumpToBtm = true, isUnderSuperRoot = false, isDummy = true
      auto * dummyRootS = new BTreeNodeT(nullptr, true, true, true, false, true);
      dummyRootS->putFirstBtm(nullptr, 0);
      srootA_.setRoot(new BTreeNodeT(dummyRootS, true, true, true, true));
      srootA_.root_->pushbackBTreeNode(dummyRootS);
    }


    /*!
     * @brief Free/delete all allocated objects.
     */
    void clearAll() {
      if (!isReady()) { // already cleared
        return;
      }
      for (uint64_t i = 0; i < idxM2S_.size() / B; ++i) {
        delete weightVecs_[i];
      }
      memutil::safefree(weightVecs_);
      { // delete separated tree
        auto * rootS = srootA_.root_->getLmBtm_DirectJump();
        while (reinterpret_cast<uintptr_t>(rootS) != BTreeNode<B>::NOTFOUND) {
          auto * next = getNextRootS(rootS);
          delete rootS;
          rootS = next;
        }
      }
      memutil::safefree(assoc_);
      idxM2S_.changeCapacity(0);
      idxS2M_.changeCapacity(0);
      memutil::safefree(parentM_);
      memutil::safefree(parentS_);
      memutil::safefree(labelM_);
      memutil::safefree(charS_);
      memutil::safefree(idxInSiblingM_);
      memutil::safefree(idxInSiblingS_);
      memutil::safefree(numChildrenS_);

      memutil::safedelete(srootM_.root_);
      memutil::safedelete(srootA_.root_);

      traCode_ = 9;
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return (srootM_.root_ != nullptr);
    }


    /*!
     * @brief Return if given 'idxM' corresponds to valid run.
     */
    bool isValidIdxM(const uint64_t idxM) const noexcept {
      return (isReady() &&
              idxM < idxM2S_.size() &&
              (idxM % B) < weightVecs_[idxM / B]->size());
    }


    /*!
     * @brief Return |T|.
     */
    size_t getSumOfWeight() const noexcept {
      assert(isReady());

      return srootM_.root_->getSumOfWeight();
    }


    /*!
     * @brief Compute num of occ of 'ch' in T.
     */
    size_t getSumOfWeight(const uint64_t ch) const noexcept {
      assert(isReady());

      const auto * retRootS = searchCharA(ch);
      if (retRootS->isDummy() || charS_[reinterpret_cast<uintptr_t>(retRootS->getLmBtm_DirectJump())] != ch) {
        return 0;
      }
      return retRootS->getSumOfWeight();
    }


    /*!
     * @brief Get length of run corresponding to 'idxM'.
     */
    uint64_t getWeightFromIdxM(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return weightVecs_[idxM / B]->read(idxM % B);
    }


    /*!
     * @brief Get character of run corresponding to 'idxM'.
     */
    uint64_t getCharFromIdxM(const uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return charS_[idxM2S_.read(idxM) / B];
    }


    /*!
     * @brief Get "idxS" from "idxM".
     */
    uint64_t idxM2S(const uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return idxM2S_.read(idxM);
    }


    /*!
     * @brief Get "idxM" from "idxS".
     */
    uint64_t idxS2M(const uint64_t idxS) const noexcept {
      return idxS2M_.read(idxS);
    }


    /*!
     * @brief Get parentS_[btmS].
     */
    auto parentS(const uint64_t btmS) const noexcept {
      return parentS_[btmS];
    }


    /*!
     * @brief Get rootM_.
     */
    auto rootM() {
      return srootM_.root_;
    }


    /*!
     * @brief Get parentS_[btmS].
     */
    auto idxInSiblingS(const uint64_t btmS) const noexcept {
      return idxInSiblingS_[btmS];
    }


    /*!
     * @brief Get associated value at "idxM".
     */
    uint64_t getAssoc(uint64_t idxM) const noexcept {
      assert(isValidIdxM(idxM));

      return assoc_[idxM];
    }


    /*!
     * @brief Set associated value at "idxM".
     */
    void setAssoc(uint64_t val, uint64_t idxM) noexcept {
      assoc_[idxM] = val;
    }


    /*!
     * @brief Get character corresponding to a node of separated tree.
     */
    uint64_t getCharFromNodeS(const BTreeNodeT * nodeS) const noexcept {
      assert(isReady());
      assert(nodeS); // nodeS should be valid node

      return charS_[reinterpret_cast<uintptr_t>(nodeS->getLmBtm_DirectJump())];
    }


    /*!
     * @brief Compute rank_{ch}[0..pos], i.e., num of ch in T[0..pos].
     */
    uint64_t rank
    (
     const uint64_t ch, //!< 64bit-char.
     uint64_t pos, //!< Pos (0base) < |T|.
     const bool calcTotalRank //!< If true, compute 'rank_{ch}[0..pos] + num of occ of characters in T smaller than ch'.
     ) const noexcept {
      assert(isReady());
      assert(pos < srootM_.root_->getSumOfWeight());

      auto idxM = searchPosM(pos); // pos is modified to relative pos
      return rank(ch, idxM, pos, calcTotalRank);
    }


    /*!
     * @brief Variant of rank function, where pos is specified by 'idxM' and 'relativePos'.
     */
    uint64_t rank
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t idxM, //!< Valid idxM.
     const uint64_t relativePos, //!< Relative pos (0base) < |T|.
     const bool calcTotalRank //!< If true, compute 'rank_{ch}[0..pos] + num of occ of characters in T smaller than ch'.
     ) const noexcept {
      assert(isValidIdxM(idxM));
      assert(relativePos < weightVecs_[idxM / B]->read(idxM % B));

      auto chNow = getCharFromIdxM(idxM);
      uint64_t ret = 0;
      uint64_t idxS;
      if (ch == chNow) {
        ret = relativePos + 1;
        idxS = idxM2S_.read(idxM);
      } else {
        const auto * retRootS = searchCharA(ch);
        if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
          return 0;
        }
        idxS = getPredIdxSFromIdxM(retRootS, ch, idxM);
      }
      const auto btmS = idxS / B;
      for (auto tmpIdxS = btmS * B; tmpIdxS < idxS + (ch != chNow); ++tmpIdxS) {
        ret += getWeightFromIdxS(tmpIdxS);
      }
      if (calcTotalRank) {
        BTreeNodeT * root;
        ret += parentS_[btmS]->calcPSum(idxInSiblingS_[btmS], root);
        return ret + root->getParent()->calcPSum(root->getIdxInSibling());
      } else {
        return ret + parentS_[btmS]->calcPSum(idxInSiblingS_[btmS]);
      }
    }


    /*!
     * @brief Compute smallest pos (0base) s.t. 'rank == rank_{ch}[0..pos]'.
     * @attention Rank is 1base.
     */
    uint64_t select
    (
     const BTreeNodeT * rootS, //!< Root of separated tree for 'ch'.
     const uint64_t rank //!< Rank > 0.
     ) const noexcept {
      assert(rank > 0);
      assert(rootS); // rootS should be valid node

      if (rank > rootS->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      }
      auto pos = rank - 1; // -1 for translating rank into 0base pos.
      const auto idxS = searchPosS(pos, rootS); // pos is modified to the relative pos
      const auto idxM = idxS2M_.read(idxS);
      const auto btmM = idxM / B;
      for (auto tmpIdxM = btmM * B; tmpIdxM < idxM; ++tmpIdxM) {
        pos += getWeightFromIdxM(tmpIdxM);
      }
      return pos + parentM_[btmM]->calcPSum(idxInSiblingM_[btmM]);
    }


    /*!
     * @brief Compute smallest pos s.t. 'rank == rank_{ch}[0..pos]'.
     * @attention Rank is 1base.
     */
    uint64_t select
    (
     const uint64_t ch, //!< character for select query.
     const uint64_t rank //!< Rank > 0.
     ) const noexcept {
      assert(rank > 0);

      const auto * retRootS = searchCharA(ch);
      if (retRootS->isDummy() || getCharFromNodeS(retRootS) != ch) {
        return BTreeNodeT::NOTFOUND;
      }
      return select(retRootS, rank);
    }


    /*!
     * @brief Compute smallest pos s.t. 'totalRank == totalRank_{ch}[0..pos]'.
     * @attention TotalRank is 1base.
     */
    uint64_t select
    (
     const uint64_t totalRank //!< TotalRank > 0.
     ) const noexcept {
      assert(totalRank > 0);

      if (totalRank > srootA_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      }
      auto pos = totalRank - 1;
      const auto * retRootS = searchPosA(pos);
      return select(retRootS, pos + 1); // +1 for 1base rank
    }


    /*!
     * @brief Output string represented by current RLE to std::ofstream.
     */
    void printString(std::ofstream & ofs) const noexcept {
      assert(isReady());

      uint64_t pos = 0;
      for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
        const size_t exponent = getWeightFromIdxM(idxM);
        char ch = getCharFromIdxM(idxM);
        for (size_t i = 0; i < exponent; ++i) {
          ofs.put(ch);
        }
      }
    }


    //// Public search functions
  public:
    /*!
     * @brief Return 'idxM' corresponding to the run containing 'pos'-th character (0base).
     * @attention 'pos' is modified to be the relative position (0base) from the beginning of the run.
     */
    uint64_t searchPosM
    (
     uint64_t & pos //!< [in,out] Give position to search (< |T|). It is modified to relative position.
     ) const noexcept {
      assert(isReady());
      assert(pos < srootM_.root_->getSumOfWeight());

      uint64_t btmM = reinterpret_cast<uintptr_t>(srootM_.root_->searchPos(pos));

      const auto * wVec = weightVecs_[btmM];
      uint8_t i = 0;
      while (pos >= wVec->read(i)) {
        pos -= wVec->read(i);
        ++i;
      }
      return btmM * B + i;
    }


    /*!
     * @brief Search root of separated tree of the largest character that is smaller or equal to 'ch'.
     */
    BTreeNodeT * searchCharA
    (
     const uint64_t ch
     ) const noexcept {
      assert(isReady());

      auto * nodeA = srootA_.root_;
      while (true) {
        const bool nowOnBorder = nodeA->isBorder();
        uint8_t lb = 0;
        uint8_t ub = nodeA->getNumChildren();
        while (lb+1 != ub) { // invariant: the answer is in [lb..ub)
          uint8_t mid = (lb + ub) / 2;
          if (ch < getCharFromNodeA(nodeA->getChildPtr(mid), nowOnBorder)) {
            ub = mid;
          } else {
            lb = mid;
          }
        }
        nodeA = nodeA->getChildPtr(lb);
        if (nowOnBorder) {
          return nodeA;
        }
      }
    }


    uint64_t searchPosS(uint64_t & pos, const BTreeNodeT * rootS) const noexcept {
      assert(isReady());
      assert(rootS); // rootS should be valid node
      assert(pos < rootS->getSumOfWeight());

      uint64_t idxS = B * reinterpret_cast<uintptr_t>(rootS->searchPos(pos));

      while (true) {
        auto weight = getWeightFromIdxS(idxS);
        if (pos >= weight) {
          pos -= weight;
          ++idxS;
        } else {
          return idxS;
        }
      }
    }


    /*!
     * Search idxS having the largest label that is smaller or equal to 'label'
     */
    uint64_t searchLabelS(const uint64_t label, const BTreeNodeT * rootS) const noexcept {
      assert(isReady());
      assert(rootS); // rootS should be valid node

      const auto * nodeS = rootS;
      while (true) {
        const bool nowOnBorder = nodeS->isBorder();
        uint8_t lb = 0;
        uint8_t ub = nodeS->getNumChildren();
        while (lb+1 != ub) {
          uint8_t mid = (lb + ub) / 2;
          if (label < getLabelFromNodeU(nodeS->getChildPtr(mid), nowOnBorder)) {
            ub = mid;
          } else {
            lb = mid;
          }
        }
        nodeS = nodeS->getChildPtr(lb);
        if (nowOnBorder) {
          break;
        }
      }
      const uint64_t idxS = B * reinterpret_cast<uintptr_t>(nodeS);
      uint8_t lb = 0;
      uint8_t ub = numChildrenS_[idxS / B];
      while (lb+1 != ub) {
        uint8_t mid = (lb + ub) / 2;
        if (label < labelM_[idxS2M_.read(idxS + mid) / B]) {
          ub = mid;
        } else {
          lb = mid;
        }
      }
      return idxS + lb;
    }


    //// Iterator like functions
  public:
    /*!
     * @brief Get previous idxM.
     */
    uint64_t getPrevIdxM
    (
     const uint64_t idxM //!< Valid idxM.
     ) const noexcept {
      assert(isValidIdxM(idxM));

      if (idxM % B) {
        return idxM - 1;
      }
      const uint64_t prevBtmM
        = reinterpret_cast<uintptr_t>(parentM_[idxM / B]->getPrevBtm(idxInSiblingM_[idxM / B]));
      if (prevBtmM != BTreeNodeT::NOTFOUND) {
        return prevBtmM * B + getNumChildrenM(prevBtmM) - 1;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get next idxM.
     */
    uint64_t getNextIdxM
    (
     const uint64_t idxM //!< Valid idxM.
     ) const noexcept {
      assert(isValidIdxM(idxM));

      if ((idxM % B) + 1 < getNumChildrenM(idxM / B)) {
        return idxM + 1;
      }
      const uint64_t nextBtmM
        = reinterpret_cast<uintptr_t>(parentM_[idxM / B]->getNextBtm_DirectJump(idxInSiblingM_[idxM / B]));
      if (nextBtmM != BTreeNodeT::NOTFOUND) {
        return nextBtmM * B;
      }
      return BTreeNodeT::NOTFOUND;
    }


    /*!
     * @brief Get first root of separated tree, which is dummy.
     */
    BTreeNodeT * getFstRootS() const noexcept {
      assert(isReady());

      return getNextRootS(srootA_.root_->getLmBtm_DirectJump());
    }


    /*!
     * @brief Get root of separated tree for previous character.
     */
    BTreeNodeT * getPrevRootS(const BTreeNodeT * node) const noexcept {
      assert(isReady());
      assert(node); // rootS should be valid node

      uint8_t idxInSib;
      do {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      } while (!node->isBorder());
      return node->getPrevBtm(idxInSib);
    }


    /*!
     * @brief Get root of separated tree for next character.
     */
    BTreeNodeT * getNextRootS(const BTreeNodeT * node) const noexcept {
      assert(isReady());
      assert(node); // rootS should be valid node

      uint8_t idxInSib;
      do {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      } while (!(node->isBorder()));
      return node->getNextBtm_DirectJump(idxInSib);
    }


    uint64_t getPrevIdxS(const size_t idxS) const noexcept {
      if (idxS % B) {
        return idxS - 1;
      }
      const uint64_t prevBtmS
        = reinterpret_cast<uintptr_t>(parentS_[idxS / B]->getPrevBtm(idxInSiblingS_[idxS / B]));
      if (prevBtmS != BTreeNodeT::NOTFOUND) {
        return prevBtmS * B + numChildrenS_[prevBtmS] - 1;
      }
      return BTreeNodeT::NOTFOUND;
    }


    uint64_t getNextIdxS(const size_t idxS) const noexcept {
      if ((idxS % B) + 1 < numChildrenS_[idxS / B]) {
        return idxS + 1;
      }
      const uint64_t nextBtmS
        = reinterpret_cast<uintptr_t>(parentS_[idxS / B]->getNextBtm_DirectJump(idxInSiblingS_[idxS / B]));
      if (nextBtmS != BTreeNodeT::NOTFOUND) {
        return nextBtmS * B;
      }
      return BTreeNodeT::NOTFOUND;
    }


    //// private getter functions (utilities)
    uint8_t getNumChildrenM(const uint64_t btmM) const noexcept {
      return weightVecs_[btmM]->size();
    }


    uint64_t getWeightFromIdxS(uint64_t idxS) const noexcept {
      const uint64_t idxM = idxS2M_.read(idxS);
      return weightVecs_[idxM / B]->read(idxM % B);
    }


    uint64_t getLabelFromNodeU(const BTreeNodeT * nodeU, const bool isChildOfBorder)  const noexcept {
      uint64_t idxS;
      if (!isChildOfBorder) {
        idxS = B * reinterpret_cast<uintptr_t>(nodeU->getLmBtm_DirectJump());
      } else {
        idxS = B * reinterpret_cast<uintptr_t>(nodeU);
      }
      uint64_t idxM = idxS2M_.read(idxS); // idxM corresponding to the left most idxS in btmS
      return labelM_[idxM / B];
    }


    uint64_t getCharFromNodeA(const BTreeNodeT * nodeA, const bool isChildOfBorder) const noexcept {
      uint64_t btmS;
      if (!isChildOfBorder) {
        btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm_DirectJump()->getLmBtm_DirectJump());
      } else {
        btmS = reinterpret_cast<uintptr_t>(nodeA->getLmBtm_DirectJump());
      }
      return charS_[btmS];
    }


    uint64_t getPrevBtmM(const uint64_t btmM) const noexcept {
      return reinterpret_cast<uintptr_t>(parentM_[btmM]->getPrevBtm(idxInSiblingM_[btmM]));
    }


    uint64_t getNextBtmM(const uint64_t btmM) const noexcept {
      return reinterpret_cast<uintptr_t>(parentM_[btmM]->getNextBtm_DirectJump(idxInSiblingM_[btmM]));
    }


    /*!
     * @brief Return root of separated tree that contains the position 'pos' (0based) in alphabetically sorted array
     */
    BTreeNodeT * searchPosA(uint64_t & pos) const noexcept {
      return srootA_.root_->searchPos(pos);
    }


    void reserveBtms(const size_t numBtms) {
      memutil::realloc_AbortOnFail(assoc_, numBtms * B);
      const uint8_t w = bits::bitSize(numBtms * B - 1);
      idxM2S_.convert(w, numBtms * B);
      idxS2M_.convert(w, numBtms * B);
      memutil::realloc_AbortOnFail(weightVecs_, numBtms);
      memutil::realloc_AbortOnFail(parentM_, numBtms);
      memutil::realloc_AbortOnFail(parentS_, numBtms);
      memutil::realloc_AbortOnFail(labelM_, numBtms);
      memutil::realloc_AbortOnFail(charS_, numBtms);
      memutil::realloc_AbortOnFail(idxInSiblingM_, numBtms);
      memutil::realloc_AbortOnFail(idxInSiblingS_, numBtms);
      memutil::realloc_AbortOnFail(numChildrenS_, numBtms);
      traCode_ = TagRelabelAlgo::getSmallestTraCode(numBtms);
    }


    void expandBtms() {
      const uint64_t newNumBtms = 2 * (idxM2S_.capacity() / B); // number of capacity of bottoms is doubled
      reserveBtms(newNumBtms);
    }


    ////
    void changeWeight(const uint64_t idxM, const int64_t change) {
      // update Leaf
      auto * wVec = weightVecs_[idxM / B];
      const uint64_t val = wVec->read(idxM % B) + change;
      const uint8_t w = wVec->getW();
      const uint8_t needW = bits::bitSize(val);
      if (needW > w) {
        wVec->convert(needW, B);
      }
      wVec->write(val, idxM % B);
      // update mixed tree
      parentM_[idxM / B]->changePSumFrom(idxInSiblingM_[idxM / B], change);
      // update separated tree AND alphabet tree (they are connected seamlessly)
      const uint64_t btmS = idxM2S_.read(idxM) / B;
      parentS_[btmS]->changePSumFrom(idxInSiblingS_[btmS], change);
    }


    void asgnLabel(const uint64_t btmM) {
      uint64_t next = getNextBtmM(btmM);
      uint64_t prev = getPrevBtmM(btmM); // assume that prev alwarys exists
      uint64_t base = (next == BTreeNodeT::NOTFOUND) ? TagRelabelAlgo::MAX_LABEL : labelM_[next];
      if (labelM_[prev] < base - 1) {
        labelM_[btmM] = (labelM_[prev] + base) / 2;
        return;
      }

      base >>= 1;
      uint64_t tmpBtmM = btmM;
      uint8_t l = 1;
      uint64_t num = 1;
      uint64_t overflowNum = 2;
      while (true) {
        while (prev != BTreeNodeT::NOTFOUND && (labelM_[prev] >> l) == base) { // expand backward
          ++num;
          tmpBtmM = prev;
          prev = getPrevBtmM(prev);
        }
        while (next != BTreeNodeT::NOTFOUND && (labelM_[next] >> l) == base){ // expand forward
          ++num;
          next = getNextBtmM(next);
        }
        if (overflowNum >= num) {
          break;
        }
        ++l;
        base >>= 1;
        overflowNum = TagRelabelAlgo::getNextOverflowNum(overflowNum, traCode_);
      }

      // relabel num labels
      uint64_t tmpLabel = base << l;
      const uint64_t interval = (UINT64_C(1) << l) / num;
      while (true) {
        labelM_[tmpBtmM] = tmpLabel;
        if (--num == 0) {
          return;
        }
        tmpLabel += interval;
        tmpBtmM = getNextBtmM(tmpBtmM);
      }
    }


    /*!
     * @brief Split btmM.
     * @post
     *   This function will do the following:
     *   - setup
     *     - parentM_[retBtmM] (by handleSplitOfBtmInBtm())
     *     - idxInSiblingM_[retBtmM] (by handleSplitOfBtmInBtm())
     *     - labelM_[retBtmM] (by asgnLabel())
     *   - resize
     *     - idxM2S_ to use range [endIdxM, endIdxM + B)
     *   - reserve
     *     - weightVecs_[retBtmM]
     *   - update
     *     - upper nodes (through handleSplitOfBtmInBtm())
     *     - labels (by asgnLabel())
     */
    uint64_t splitBtmM(const uint8_t width, const uint64_t btmM, const uint64_t weight) {
      const uint64_t endIdxM = idxM2S_.size();
      const uint64_t retBtmM = endIdxM / B;
      if (!(idxM2S_.resizeWithoutReserve(endIdxM + B))) {
        expandBtms();
      }
      idxM2S_.resize(endIdxM + B);
      // reserve
      weightVecs_[retBtmM] = new WBitsVec(width, B);
      // setup and update
      handleSplitOfBtmInBtm(btmM, retBtmM, weight, parentM_, idxInSiblingM_);
      asgnLabel(retBtmM);
      return retBtmM;
    }


    /*!
     * @brief Split btmS.
     * @post
     *   This function will do the following:
     *   - setup
     *     - parentS_[retBtmS] (by handleSplitOfBtmInBtm())
     *     - idxInSiblingS_[retBtmS] (by handleSplitOfBtmInBtm())
     *     - charS_[retBtmS]
     *   - resize
     *     - idxS2M_ to use range [endIdxS, endIdxS + B)
     *   - update
     *     - upper nodes (through handleSplitOfBtmInBtm())
     */
    uint64_t splitBtmS(const uint64_t btmS, const uint64_t weight) {
      const uint64_t endIdxS = idxS2M_.size();
      const uint64_t retBtmS = endIdxS / B;
      if (!(idxS2M_.resizeWithoutReserve(endIdxS + B))) {
        expandBtms();
      }
      idxS2M_.resize(endIdxS + B);
      // setup and update
      handleSplitOfBtmInBtm(btmS, retBtmS, weight, parentS_, idxInSiblingS_);
      charS_[retBtmS] = charS_[btmS];
      return retBtmS;
    }


    uint64_t setupNewSTree(BTreeNodeT * predNode, const uint64_t ch) {
      const uint64_t endIdxS = idxS2M_.size();
      const uint64_t btmS = endIdxS / B;
      if (!(idxS2M_.resizeWithoutReserve(endIdxS + B))) {
        expandBtms();
      }
      idxS2M_.resize(endIdxS + B);
    
      auto * newRootS = new BTreeNodeT(reinterpret_cast<BTreeNodeT *>(btmS), true, true, true, false);
      parentS_[btmS] = newRootS;
      idxInSiblingS_[btmS] = 0;
      charS_[btmS] = ch;
      numChildrenS_[btmS] = 1; // only dummy idxS exists
      idxS2M_.write(0, btmS * B); // link to dummy idxM of weight 0

      newRootS->pushbackBtm(reinterpret_cast<BTreeNodeT *>(btmS), 0);
      const auto idxInSib = predNode->getIdxInSibling();
      auto * parent = predNode->getParent();
      parent->handleSplitOfChild(newRootS, idxInSib);
      return endIdxS;
    }


    void mvIdxFwd(WBitsVec & wba, uint64_t srcIdx, uint64_t tgtIdx, uint64_t num, WBitsVec & wbaOther) {
      for (uint64_t i = num; i > 0; --i) {
        const uint64_t idxOther = wba.read(srcIdx + i - 1);
        wba.write(idxOther, tgtIdx + i - 1);
        wbaOther.write(tgtIdx + i - 1, idxOther);
      }
    }


    uint64_t makeSpaceAfterIdxM(const uint64_t idxM) {
      const uint8_t remIdxM = idxM % B;
      const uint64_t btmM = idxM / B;
      WBitsVec * wVec0 = weightVecs_[btmM];
      const uint8_t oriNum = wVec0->size();
      if (oriNum < B) {
        wVec0->resize(oriNum + 1);
        const uint8_t mvNum = oriNum - remIdxM - 1;
        if (mvNum) {
          mvWBA_SameW(wVec0->getItrAt(remIdxM + 1), wVec0->getItrAt(remIdxM + 2), mvNum);
          mvIdxFwd(idxM2S_, idxM + 1, idxM + 2, mvNum, idxS2M_);
          bits::cpBytes(assoc_ + idxM + 1, assoc_ + idxM + 2, mvNum * (sizeof(assoc_[0])));
        }
        wVec0->write(0, remIdxM + 1);
        return idxM + 1;
      }
      // split
      uint64_t sum = 0;
      for (uint8_t i = B/2; i < B; ++i) {
        sum += wVec0->read(i);
      }
      const auto newBtmM = splitBtmM(wVec0->getW(), btmM, sum);
      WBitsVec * wVec1 = weightVecs_[newBtmM];
      mvWBA_SameW(wVec0->getItrAt(B/2), wVec1->getItrAt(0), B/2);
      wVec0->resize(B/2);
      wVec1->resize(B/2);
      mvIdxFwd(idxM2S_, btmM*B + B/2, newBtmM*B, B/2, idxS2M_);
      bits::cpBytes(assoc_ + btmM*B + B/2, assoc_ + newBtmM*B, B/2 * (sizeof(assoc_[0])));
      if (remIdxM < B/2) {
        return makeSpaceAfterIdxM(idxM);
      } else {
        return makeSpaceAfterIdxM(newBtmM * B + remIdxM - B/2);
      }
    }


    uint64_t makeSpaceAfterIdxS(const uint64_t idxS) {
      const uint8_t remIdxS = idxS % B;
      const uint64_t btmS = idxS / B;
      const uint8_t oriNum = numChildrenS_[btmS];
      if (oriNum < B) {
        numChildrenS_[btmS] = oriNum + 1;
        const uint8_t mvNum = oriNum - remIdxS - 1;
        if (mvNum) {
          mvIdxFwd(idxS2M_, idxS + 1, idxS + 2, mvNum, idxM2S_);
        }
        return idxS + 1;
      }
      uint64_t sum = 0;
      for (uint8_t i = B/2; i < B; ++i) {
        sum += getWeightFromIdxS(btmS * B + i);
      }
      const auto newBtmS = splitBtmS(btmS, sum);
      numChildrenS_[btmS] = B/2;
      numChildrenS_[newBtmS] = B/2;
      mvIdxFwd(idxS2M_, btmS*B + B/2, newBtmS*B, B/2, idxM2S_);
      if (remIdxS < B/2) {
        return makeSpaceAfterIdxS(idxS);
      } else {
        return makeSpaceAfterIdxS(newBtmS * B + remIdxS - B/2);
      }
    }


    void handleSplitOfBtmInBtm
    (
     const uint64_t btm,
     const uint64_t newBtm,
     const uint64_t weight,
     BTreeNodeT ** parentArray,
     uint8_t * idxInSibArray
     ) {
      auto * uNode = parentArray[btm];
      const auto idxInSib = idxInSibArray[btm];
      const auto oriNum = uNode->getNumChildren();
      const uint8_t numToL = uNode->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(newBtm), weight, idxInSib);
      if (numToL == 0) {
        for (uint8_t i = idxInSib + 1; i < uNode->getNumChildren(); ++i) {
          const uint64_t tmp = reinterpret_cast<uintptr_t>(uNode->getChildPtr(i));
          parentArray[tmp] = uNode;
          idxInSibArray[tmp] = i;
        }
        if (oriNum == B) {
          auto * nextNode = uNode->getNextSib();
          for (uint8_t i = 0; i < nextNode->getNumChildren(); ++i) {
            const uint64_t tmp = reinterpret_cast<uintptr_t>(nextNode->getChildPtr(i));
            parentArray[tmp] = nextNode;
            idxInSibArray[tmp] = i;
          }
        }
      } else {
        for (uint8_t i = 0; i < uNode->getNumChildren(); ++i) {
          const uint64_t tmp = reinterpret_cast<uintptr_t>(uNode->getChildPtr(i));
          parentArray[tmp] = uNode;
          idxInSibArray[tmp] = i;
        }
        auto * prevNode = uNode->getPrevSib();
        const uint8_t lnum = prevNode->getNumChildren();
        for (uint8_t i = lnum - (numToL + (idxInSib < numToL)); i < lnum; ++i) {
          const uint64_t tmp = reinterpret_cast<uintptr_t>(prevNode->getChildPtr(i));
          parentArray[tmp] = prevNode;
          idxInSibArray[tmp] = i;
        }
      }
    }


    uint64_t getPredIdxSFromIdxM(const BTreeNodeT * rootS, const uint64_t ch, const uint64_t idxM) const noexcept {
      const uint64_t btmM = idxM / B;
      if (btmM) { // If btmM is not 0 (0 means btmM is the first btm in the mixed tree).
        uint64_t i = idxM - 1;
        for ( ; i >= btmM * B && getCharFromIdxM(i) != ch; --i) {}
        if (i >= btmM * B) {
          return idxM2S_.read(i);
        } else {
          return searchLabelS(labelM_[btmM] - 1, rootS); // -1 is needed.
        }
      } else { // btmM == 0: dummy idx (== 0) should be ignored.
        uint64_t i = idxM - 1;
        for ( ; i > 0 && getCharFromIdxM(i) != ch; --i) {}
        if (i > 0) {
          return idxM2S_.read(i);
        } else {
          return B * reinterpret_cast<uintptr_t>(rootS->getLmBtm_DirectJump());
        }
      }
    }


    /*!
     * @brief Insert new run of character 'ch' and length 'weight' after 'idxM'.
     * @return IdxM of the inserted run.
     */
    uint64_t insertNewRunAfter(const uint64_t ch, const uint64_t weight, const uint64_t idxM) {
      const auto newIdxM = makeSpaceAfterIdxM(idxM);
      auto * retRootS = searchCharA(ch);
      uint64_t idxS;
      if (retRootS->isDummy() || charS_[reinterpret_cast<uintptr_t>(retRootS->getLmBtm_DirectJump())] != ch) {
        idxS = setupNewSTree(retRootS, ch);
      } else {
        idxS = getPredIdxSFromIdxM(retRootS, ch, newIdxM);
      }
      const auto newIdxS = makeSpaceAfterIdxS(idxS);
      idxM2S_.write(newIdxS, newIdxM);
      idxS2M_.write(newIdxM, newIdxS);
      changeWeight(newIdxM, weight);
      return newIdxM;
    }


  public:
    /*!
     * @brief Pushback a run, merging into the last run if possible.
     */
    uint64_t pushbackRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     uint64_t & pos //!< [out] It is set to relative position of a run.
     ) {
      const uint64_t btmM = reinterpret_cast<uintptr_t>(srootM_.root_->getRmBtm());
      const auto idxM = btmM * B + getNumChildrenM(btmM) - 1;
      if (getCharFromIdxM(idxM) != ch) {
        pos = 0;
        return insertNewRunAfter(ch, weight, idxM);
      } else { // Merge into the last run
        pos = getWeightFromIdxM(idxM);
        changeWeight(idxM, weight);
        return idxM;
      }
    }


    /*!
     * @brief Pushback a run without merge.
     */
    uint64_t pushbackRunWithoutMerge
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight //!< Weight (exponent) of new run.
     ) {
      const uint64_t btmM = reinterpret_cast<uintptr_t>(srootM_.root_->getRmBtm());
      return insertNewRunAfter(ch, weight, btmM * B + getNumChildrenM(btmM) - 1);
    }


    /*!
     * @brief Insert run of 'ch^{weight}' at 'pos', merging into adjacent runs if possible.
     */
    uint64_t insertRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     uint64_t & pos //!< [in,out] 0base position where inserted run will start. It is modified to relative position in a run.
     ) {
      if (pos > srootM_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      } else if (pos == srootM_.root_->getSumOfWeight()) {
        return pushbackRun(ch, weight, pos);
      }
      auto idxM = searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      auto chNow = getCharFromIdxM(idxM);
      if (ch == chNow) {
        changeWeight(idxM, weight);
      } else if (pos == 0) {
        idxM = getPrevIdxM(idxM); // Move to previous idxM.
        if (idxM > 0 && ch == getCharFromIdxM(idxM)) { // Check if 'ch' can be merged with the previous run.
          pos = getWeightFromIdxM(idxM);
          changeWeight(idxM, weight);
        } else {
          idxM = insertNewRunAfter(ch, weight, idxM);
        }
      } else { // Current run is split with fstHalf of weight 'pos'.
        const auto weightSndHalf = getWeightFromIdxM(idxM) - pos;
        pos = 0;
        changeWeight(idxM, -1 * weightSndHalf);
        idxM = insertNewRunAfter(ch, weight, idxM);
        idxM = insertNewRunAfter(chNow, weightSndHalf, idxM);
        idxM = getPrevIdxM(idxM);
      }
      return idxM;
    }


    /*!
     * @brief Variant of DynRLE::insertRun for rvalue pos.
     */
    uint64_t insertRun
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     uint64_t && pos //!< 0base position where inserted run will start.
     ) {
      auto tmp = pos;
      return insertRun(ch, weight, tmp);
    }


    /*!
     * @brief Insert run of 'ch^{weight}' at 'pos' without merge.
     */
    uint64_t insertRunWithoutMerge
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     uint64_t & pos //!< [in,out] 0base position where inserted run will start.
     ) {
      if (pos > srootM_.root_->getSumOfWeight()) {
        return BTreeNodeT::NOTFOUND;
      } else if (pos == srootM_.root_->getSumOfWeight()) {
        pos = 0;
        return pushbackRunWithoutMerge(ch, weight);
      }
      auto idxM = searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
      if (pos != 0) { // Current run is split with fstHalf of weight 'pos'.
        auto chNow = getCharFromIdxM(idxM);
        const auto weightSndHalf = getWeightFromIdxM(idxM) - pos;
        changeWeight(idxM, -1 * weightSndHalf);
        idxM = insertNewRunAfter(chNow, weightSndHalf, idxM);
      }
      idxM = getPrevIdxM(idxM); // Move to previous idxM.
      idxM = insertNewRunAfter(ch, weight, idxM);
      pos = 0;
      return idxM;
    }


    /*!
     * @brief Variant of DynRLE::insertRunWithoutMerge for rvalue pos.
     */
    uint64_t insertRunWithoutMerge
    (
     const uint64_t ch, //!< 64bit-char.
     const uint64_t weight, //!< Weight (exponent) of new run.
     uint64_t && pos //!< 0base position where inserted run will start.
     ) {
      auto tmp = pos;
      return insertRunWithoutMerge(ch, weight, tmp);
    }


    //// statistics
  public:
    size_t calcMemBytesAssoc() const noexcept {
      return idxM2S_.capacity() * sizeof(assoc_[0]);
    }


    size_t calcMemBytesMTree() const noexcept {
      return srootM_.root_->calcMemBytes();
    }


    size_t calcMemBytesATree() const noexcept {
      return srootA_.root_->calcMemBytes();
    }


    size_t calcMemBytesSTree() const noexcept {
      size_t size = 0;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        size += rootS->calcMemBytes();
      }
      return size;
    }


    size_t calcMemBytesWeightVecs() const noexcept {
      size_t size = 0;
      for (uint64_t i = 0; i < idxM2S_.size() / B; ++i) {
        size += weightVecs_[i]->calcMemBytes();
      }
      return size;
    }


    size_t calcMemBytesIdxConvertVecs() const noexcept {
      size_t size = 0;
      size += idxM2S_.calcMemBytes();
      size += idxS2M_.calcMemBytes();
      return size;
    }


    size_t calcMemBytesBtmArrays() const noexcept {
      return (idxM2S_.capacity() / B) * (sizeof(parentM_[0]) + sizeof(parentS_[0]) +
                                         sizeof(labelM_[0]) + sizeof(charS_[0]) +
                                         sizeof(idxInSiblingM_[0]) + sizeof(idxInSiblingS_[0]) +
                                         sizeof(numChildrenS_[0]) + sizeof(weightVecs_[0]));
    }


    size_t calcMemBytes() const noexcept {
      size_t size = sizeof(*this);
      size += calcMemBytesAssoc();
      size += calcMemBytesMTree();
      size += calcMemBytesATree();
      size += calcMemBytesSTree();
      size += calcMemBytesWeightVecs();
      size += calcMemBytesIdxConvertVecs();
      size -= sizeof(idxM2S_); // minus double counted part
      size -= sizeof(idxS2M_); // minus double counted part
      size += calcMemBytesBtmArrays();
      return size;
    }


    size_t calcNumUsedSTree() const noexcept {
      size_t numUsed = 0;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        numUsed += rootS->calcNumUsed();
      }
      return numUsed;
    }


    size_t calcNumSlotsSTree() const noexcept {
      size_t numSlots = 0;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        numSlots += rootS->calcNumSlots();
      }
      return numSlots;
    }


    size_t calcNumRuns() const noexcept {
      size_t numRuns = 0;
      for (size_t i = 0; i < idxM2S_.size() / B; ++i) {
        numRuns += getNumChildrenM(i);
      }
      return numRuns - 1; // -1 due to the first dummy
    }


    size_t calcNumAlph() const noexcept {
      size_t numAlph = 0;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        ++numAlph;
      }
      return numAlph;
    }


    void printStatictics(std::ostream & os) const noexcept {
      const size_t totalLen = getSumOfWeight();
      const size_t numRuns = calcNumRuns();
      os << "TotalLen = " << totalLen << ", #Runs = " << numRuns << ", Alphabet Size = " << calcNumAlph() << ", BTree arity param B = " << static_cast<int>(B) << std::endl;
      os << "Total: " << calcMemBytes() << " bytes" << std::endl;
      os << "Assoc: " << calcMemBytesAssoc() << " bytes" << std::endl;
      os << "MTree: " << calcMemBytesMTree() << " bytes, OccuRate = " << ((srootM_.root_->calcNumSlots()) ? 100.0 * srootM_.root_->calcNumUsed() / srootM_.root_->calcNumSlots() : 0)
         << " (= 100*" << srootM_.root_->calcNumUsed() << "/" << srootM_.root_->calcNumSlots() << ")" << std::endl;
      os << "ATree: " << calcMemBytesATree() << " bytes, OccuRate = " << ((srootA_.root_->calcNumSlots()) ? 100.0 * srootA_.root_->calcNumUsed() / srootA_.root_->calcNumSlots() : 0)
         << " (= 100*" << srootA_.root_->calcNumUsed() << "/" << srootA_.root_->calcNumSlots() << ")" << std::endl;
      os << "STree: " << calcMemBytesSTree() << " bytes, OccuRate = " << ((calcNumSlotsSTree()) ? 100.0 * calcNumUsedSTree() / calcNumSlotsSTree() : 0)
         << " (= 100*" << calcNumUsedSTree() << "/" << calcNumSlotsSTree() << ")" << std::endl;
      os << "IdxConvertVecs: " << calcMemBytesIdxConvertVecs() << " bytes ~ "
         << "(2*" << static_cast<int>(idxM2S_.getW()) << "(bitwidth)*" << idxM2S_.capacity() << "(capacity each))/8, "
         << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * 2 * numRuns / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
         << " (= 100*2*" << numRuns << "/" << (idxM2S_.capacity() + idxS2M_.capacity()) << ")" << std::endl;
      os << "WeightVecs: " << calcMemBytesWeightVecs() << " bytes" << std::endl;
      os << "BtmArrays: " << calcMemBytesBtmArrays() << " bytes, "
         << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * (idxM2S_.size() + idxS2M_.size()) / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
         << " (= 100*" << (idxM2S_.size() + idxS2M_.size())/B << "/" << (idxM2S_.capacity() + idxS2M_.capacity())/B << "), "
         << "OccuRate (btmM) = " << ((idxM2S_.capacity()) ? 100.0 * idxM2S_.size() / idxM2S_.capacity() : 0)
         << " (= 100*" << idxM2S_.size()/B << "/" << idxM2S_.capacity()/B << "), "
         << "OccuRate (btmS) = " << ((idxS2M_.capacity()) ? 100.0 * idxS2M_.size() / idxS2M_.capacity() : 0)
         << " (= 100*" << idxS2M_.size()/B << "/" << idxS2M_.capacity()/B << ")" << std::endl;
    }


    void printDebugInfo(std::ostream & os) const noexcept {
      { // check links of idxM2S and idxS2M
        const uint64_t numBtmM = idxM2S_.size() / B;
        for (uint64_t i = 0; i < numBtmM; ++i) {
          for (uint64_t j = 0; j < getNumChildrenM(i); ++j) {
            if (j < getNumChildrenM(i) && B*i+j != idxS2M_.read(idxM2S_.read(B*i+j))) {
              os << "error!! links of idxM2S and idxS2M" << std::endl; // WARNING, links are not maintained correctly
            }
          }
        }
      }


      { // check links of parent-child for M
        const uint64_t numBtmM = idxM2S_.size() / B;
        for (uint64_t i = 0; i < numBtmM; ++i) {
          uint8_t idx = idxInSiblingM_[i];
          auto node = parentM_[i];
          bool islmbtm = (idx == 0);
          if (reinterpret_cast<uintptr_t>(node->getChildPtr(idx)) != i) {
            os << "error!! " << "parent-child for btmM = " << i << std::endl;
          }
          if (islmbtm && reinterpret_cast<uintptr_t>(node->getLmJumpNode()) != i) {
            os << "error!! lmJumNode for btmM = " << i << std::endl;
          }
          while (!(node->isRoot())) {
            idx = node->getIdxInSibling();
            islmbtm &= (idx == 0);
            if (node->getParent()->getChildPtr(idx) != node) {
              os << "error!! " << "parent-child for child node = " << node << std::endl;
            }
            if (islmbtm && reinterpret_cast<uintptr_t>(node->getLmJumpNode()) != i) {
              os << "error!! lmJumNode for btmM = " << i << std::endl;
            }
            node = node->getParent();
          }
        }
      }


      { // check links of parent-child for S
        const uint64_t numBtmS = idxS2M_.size() / B;
        for (uint64_t i = 0; i < numBtmS; ++i) {
          uint8_t idx = idxInSiblingS_[i];
          auto node = parentS_[i];
          bool islmbtm = (idx == 0);
          if (reinterpret_cast<uintptr_t>(node->getChildPtr(idx)) != i) {
            os << "error!! " << "parent-child for btmS = " << i << std::endl;
          }
          if (islmbtm && reinterpret_cast<uintptr_t>(node->getLmJumpNode()) != i) {
            os << "error!! lmJumpNode for btmS = " << i << std::endl;
          }
          while (!(node->isRoot())) {
            idx = node->getIdxInSibling();
            islmbtm &= (idx == 0);
            if (node->getParent()->getChildPtr(idx) != node) {
              os << "error!! " << "parent-child for child node = " << node << std::endl;
            }
            if (islmbtm && reinterpret_cast<uintptr_t>(node->getLmJumpNode()) != i) {
              os << "error!! lmJumNode for btmM = " << i << std::endl;
            }
            node = node->getParent();
          }
        }
      }

      { // check correctness of runs
        uint64_t c = UINT64_MAX;
        std::cout << "check runs:" << std::endl;
        // std::cout << srootM_.root_ << " " << srootA_.root_ << std::endl;
        uint64_t pos = 0;
        uint64_t len = 0;
        for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
          ++pos;
          len += getWeightFromIdxM(idxM);
          if (getWeightFromIdxM(idxM) == 0) {
            std::cout << "detected 0 length run: " << idxM << ", " << pos << std::endl;
          }
          if (c == getCharFromIdxM(idxM)) {
            auto idxM0 = getPrevIdxM(idxM);
            std::cout << "detected consecutive runs having the same char: " 
                      << idxM << ", " << pos << ", (" << c << ", " << getWeightFromIdxM(idxM0) << ")" << ", (" << c << ", " << getWeightFromIdxM(idxM) << ")" << std::endl;
          }
          c = getCharFromIdxM(idxM);
        }
        std::cout << "run: " << pos << ", len: " << len << std::endl;
      }

      {
        uint64_t pos = 0;
        for (auto idxM = searchPosM(pos); idxM != BTreeNodeT::NOTFOUND; idxM = getNextIdxM(idxM)) {
          os << "(" << idxM << ":" << getCharFromIdxM(idxM) << "^" << getWeightFromIdxM(idxM) << ", " << getAssoc(idxM) << ") ";
        }
        os << std::endl;
      }

      // {
      //   const uint64_t numBtmM = idxM2S_.size() / B;
      //   os << "information on M" << std::endl;
      //   for (uint64_t i = 0; i < numBtmM; ++i) {
      //     const auto nextBtmM = getNextBtmM(i);
      //     os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)getNumChildrenM(i) << " lbl=" 
      //        << labelM_[i] << " par=" << parentM_[i] << " sib=" << (int)idxInSiblingM_[i] << ") "
      //        << "=> " << nextBtmM * B << std::endl;
      //     for (uint64_t j = 0; j < getNumChildrenM(i); ++j) {
      //       if (j < getNumChildrenM(i) && B*i+j != idxS2M_.read(idxM2S_.read(B*i+j))) {
      //         os << "!!"; // WARNING, links are not maintained correctly
      //       }
      //       os << idxM2S_.read(B*i+j) << "(" << getWeightFromIdxM(B*i+j) << ")  ";
      //     }
      //     os << std::endl;
      //   }
      // }

      // {
      //   const uint64_t numBtmS = idxS2M_.size() / B;
      //   os << "information on S" << std::endl;
      //   for (uint64_t i = 0; i < numBtmS; ++i) {
      //     const auto nextIdxS = getNextIdxS(i*B + numChildrenS_[i] - 1);
      //     os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)numChildrenS_[i] << " ch=" << charS_[i] << " par=" 
      //        << parentS_[i] << " sib=" << (int)idxInSiblingS_[i] << ") "
      //        << "=> " << nextIdxS << std::endl;
      //     for (uint64_t j = 0; j < B; ++j) {
      //       os << idxS2M_.read(B*i+j) << "  ";
      //     }
      //     os << std::endl;
      //   }
      // }

      os << "Alphabet: " << std::endl;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNodeT::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        const uint64_t btmS = reinterpret_cast<uintptr_t>(rootS->getLmBtm_DirectJump());
        os << "(" << charS_[btmS] << ", " << rootS->getSumOfWeight() << ") ";
      }
      os << std::endl;
    }
  };






  //////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
#include <ostream>
#include <fstream>

  using bwtintvl = std::pair<uint64_t, uint64_t>;
  using bwttracker = std::tuple<uint64_t, uint64_t, uint64_t>;

  /*!
   * @brief Online Run-length encoded Burrows–Wheeler transform (RLBWT).
   * @note
   *   ::OnlineRlbwt wraps ::DynRLE to use it for representing dynamic RLE of BWT.
   *   In contrast to ::DynRle, OnlineRLWT has a vertial end marker (em_) at emPos_.
   */
  template <class DynRle, uint8_t B = 64>
  class OnlineRlbwt
  {
  public:
    using BTreeNodeT = BTreeNode<B>;


  private:
    DynRle drle_;
    uint64_t emPos_; //!< Current position (0base) of end marker.
    uint64_t em_; //!< End marker that should not appear in the input text.
    uint64_t succSamplePos_; // Tracking txt-pos (0base) for bwt-position next to current emPos_.


  public:
    OnlineRlbwt
    (
     const size_t initNumBtms, //!< Give initial size of DynRle to reserve.
     uint64_t em = UINT64_MAX //!< Give end marker (default UINT64_MAX).
     ) :
      drle_(initNumBtms),
      emPos_(0),
      em_(em),
      succSamplePos_(0)
    {}


    /*!
     * @brief Get end marker.
     */
    uint64_t getEm() const noexcept {
      return em_;
    }


    /*!
     * @brief Get current position of end marker.
     */
    uint64_t getEmPos() const noexcept {
      return emPos_;
    }


    /*!
     * @brief Get current "succSamplePos_".
     */
    uint64_t getSuccSamplePos() const noexcept {
      return succSamplePos_;
    }

  
    /*!
     * @brief Get associated value at "idxM".
     */
    uint64_t getAssoc(uint64_t idxM) const noexcept {
      return drle_.getAssoc(idxM);
    }


    /*!
     * @brief Set associated value at "idxM".
     */
    void setAssoc(uint64_t val, uint64_t idxM) noexcept {
      drle_.setAssoc(val, idxM);
    }


    /*!
     * @brief Extend RLBWT by appending one.
     */
    void extend
    (
     const uint64_t ch, //!< 64bit-char to append.
     const uint64_t txtPos //!< Txt-position of "ch" (0base).
     ) {
      uint64_t idxM;
      auto pos = emPos_;
      if (pos == drle_.getSumOfWeight()) {
        idxM = drle_.pushbackRun(ch, 1, pos);
      } else {
        idxM = drle_.searchPosM(pos); // 'pos' is modified to be the relative pos in the run of 'idxM'.
        auto chNow = drle_.getCharFromIdxM(idxM);
        if (ch == chNow) {
          drle_.changeWeight(idxM, 1);
        } else if (pos == 0) {
          idxM = drle_.getPrevIdxM(idxM); // Move to previous idxM.
          if (idxM > 0 && ch == drle_.getCharFromIdxM(idxM)) { // Check if 'ch' can be merged with the previous run.
            pos = drle_.getWeightFromIdxM(idxM);
            drle_.changeWeight(idxM, 1);
          } else {
            idxM = drle_.insertNewRunAfter(ch, 1, idxM);
          }
        } else { // Current run is split with fstHalf of weight 'pos'.
          const auto weightSndHalf = drle_.getWeightFromIdxM(idxM) - pos;
          pos = 0;
          drle_.changeWeight(idxM, -1 * weightSndHalf);
          idxM = drle_.insertNewRunAfter(ch, 1, idxM);
          idxM = drle_.insertNewRunAfter(chNow, weightSndHalf, idxM);
          drle_.setAssoc(succSamplePos_, idxM);
          idxM = drle_.getPrevIdxM(idxM);
        }
      }
      if (pos == 0) {
        drle_.setAssoc(txtPos, idxM);
      }

      if (pos + 1 != drle_.getWeightFromIdxM(idxM)) {
        ++succSamplePos_;
      } else {
        uint64_t idxS = drle_.getNextIdxS(drle_.idxM2S(idxM));
        if (idxS != DynRle::BTreeNodeT::NOTFOUND) { // Succcessor with "ch" was found.
          succSamplePos_ = drle_.getAssoc(drle_.idxS2M(idxS)) + 1;
        } else { // Succcessor with "ch" was NOT found.
          /* Take the smallest character larger "next_ch" than "ch".
             If such character exists, set "succSamplePos_" to the
             sampled position for the first BWT-run of "next_ch" minus one. */
          const auto * retRootS = drle_.searchCharA(ch);
          const auto nextRootS = drle_.getNextRootS(retRootS);
          if (reinterpret_cast<uintptr_t>(nextRootS) != DynRle::BTreeNodeT::NOTFOUND) {
            idxS = reinterpret_cast<uintptr_t>(nextRootS->getLmBtm_DirectJump()) * B + 1;
            succSamplePos_ = drle_.getAssoc(drle_.idxS2M(idxS)) + 1;
          }
        }
      }

      emPos_ = drle_.rank(ch, idxM, pos, true);
    }


    /*!
     * @brief Access to the current RLBWT by [] operator.
     */
    uint64_t operator[]
    (
     uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEm()].
     ) const noexcept {
      assert(pos < getLenWithEm());

      if (pos == emPos_) {
        return em_;
      } else if (pos > emPos_) {
        --pos;
      }
      uint64_t idxM = drle_.searchPosM(pos);
      return drle_.getCharFromIdxM(idxM);
    }


    /*!
     * @brief Return current length including end marker.
     */
    uint64_t getLenWithEm() const noexcept {
      return drle_.getSumOfWeight() + 1; // +1 for end marker, which is not in drle_.
    }


    /*!
     * @brief Return 'rank of ch at pos' + 'num of total occ of characters smaller than ch'.
     */
    uint64_t totalRank
    (
     uint64_t ch,
     uint64_t pos //!< in [0, OnlineRLBWT::getLenWithEm()].
     ) const noexcept {
      assert(pos < getLenWithEm());

      if (pos > emPos_) {
        --pos;
      }
      return drle_.rank(ch, pos, true);
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W
     * @note Intervals are [left, right) : right bound is excluded
     */
    bool lfMap
    (
     bwttracker & tracker,
     const uint64_t ch
     ) const noexcept {
      assert(ch != getEm());
      assert(std::get<0>(tracker) <= getLenWithEm() && std::get<1>(tracker) <= getLenWithEm());
      assert(std::get<0>(tracker) < std::get<1>(tracker));

      // If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch) {
        return false;
      }

      uint64_t r_in_drle = std::get<1>(tracker) - (std::get<1>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      /* +1 because in F we are not taking into account the end-marker,
         which is in position 0 but not explicitly stored in F. */
      r_in_drle = drle_.rank(ch, r_in_drle - 1, true) + 1;
      if (r_in_drle <= 1) {
        return false;
      }

      uint64_t l_in_drle = std::get<0>(tracker) - (std::get<0>(tracker) > emPos_); // Taking (implicit) end-marker into account.
      const uint64_t idxM = drle_.searchPosM(l_in_drle); // l_in_drle is modified to relative pos.
      /* Replicate variant of rank function, where pos is specified
         by 'idxM' and 'relativePos' with several modifications. */
      const auto chNow = drle_.getCharFromIdxM(idxM);
      uint64_t idxS;
      {
        if (ch == chNow) {
          idxS = drle_.idxM2S(idxM);
        } else {
          l_in_drle = 0;
          idxS = drle_.getPredIdxSFromIdxM(retRootS, ch, idxM);
        }
        const auto btmS = idxS / B;
        for (auto tmpIdxS = btmS * B; tmpIdxS < idxS + (ch != chNow); ++tmpIdxS) {
          l_in_drle += drle_.getWeightFromIdxS(tmpIdxS);
        }
        BTreeNodeT * rootS;
        l_in_drle += drle_.parentS(btmS)->calcPSum(drle_.idxInSiblingS(btmS), rootS);
        l_in_drle += rootS->getParent()->calcPSum(rootS->getIdxInSibling());
      }
      ++l_in_drle; // +1 for (implicit) end-marker in F.
    
      if (l_in_drle >= r_in_drle) {
        return false;
      }

      // Update tracker.
      std::get<0>(tracker) = l_in_drle;
      std::get<1>(tracker) = r_in_drle;
      if (l_in_drle == emPos_) {
        std::get<2>(tracker) = drle_.getSumOfWeight();
      } else if (ch == chNow) {
        std::get<2>(tracker) += 1;
      } else {
        std::get<2>(tracker) = drle_.getAssoc(drle_.idxS2M(drle_.getNextIdxS(idxS))) + 1;
      }

      return true;
    }


    /*!
     * @brief Compute bwt-interval for cW from bwt-interval for W
     * @note Intervals are [left, right) : right bound is excluded
     */
    bwtintvl lfMap
    (
     bwtintvl intvl,
     uint64_t ch
     ) const noexcept {
      assert(ch != getEm());
      assert(intvl.first <= getLenWithEm() && intvl.second <= getLenWithEm());

      // If "ch" is not in the alphabet or empty interval, return empty interval.
      const auto * retRootS = drle_.searchCharA(ch);
      if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch || intvl.first >= intvl.second) {
        return {0, 0};
      }

      uint64_t l = intvl.first - (intvl.first > emPos_);
      uint64_t r = intvl.second - (intvl.second > emPos_);
      const uint64_t idxM = drle_.searchPosM(l); // l is modified to relative pos.
      // +1 because in F.select(0, ch) we are not taking into account the end-marker,
      // which is in position 0 but not explicitly stored in F.
      return {
        drle_.rank(ch, idxM, l, true) - (drle_.getCharFromIdxM(idxM) == ch) + 1,
          drle_.rank(ch, r-1, true) + 1
          };
    }


    /*!
     * @brief LF map.
     */
    uint64_t lfMap(uint64_t i){
      assert(i < getLenWithEm());

      if (i > emPos_) {
        --i;
      }
      const uint64_t idxM = drle_.searchPosM(i);
      const unsigned char ch = drle_.getCharFromIdxM(idxM);
      return drle_.rank(ch, idxM, i, true);
    }


    /*!
     * @brief Print statistics of ::DynRLE (not of ::OnlineRLBWT).
     */
    void printStatictics
    (
     std::ostream & os //!< std::ostream (e.g., std::cout).
     ) const noexcept {
      drle_.printStatictics(os);
    }


    void printDebugInfo(std::ostream & os) const noexcept {
      std::cout <<  "emPos_ = " << emPos_ << ", em_ = " << em_ << ", succSamplePos_ = " << succSamplePos_ << std::endl;
      drle_.printDebugInfo(os);
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes() const noexcept {
      return sizeof(*this) + drle_.calcMemBytes();
    }


    /*!
     * @brief Output original text to std::ofstream.
     */
    void invert
    (
     std::ofstream & ofs
     ) const noexcept {
      uint64_t pos = 0;
      for (uint64_t i = 0; i < this->getLenWithEm() - 1; ++i) {
        if (pos > emPos_) {
          --pos;
        }
        const uint64_t idxM = drle_.searchPosM(pos);
        const unsigned char ch = drle_.getCharFromIdxM(idxM);
        ofs.put(ch);
        pos = drle_.rank(ch, idxM, pos, true);
      }
    }
  };

} // namespace itmmti

#endif
