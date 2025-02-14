/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * Modified by [Xinwu Yu] [2025]
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file OnlineRlbwt.hpp
 * @brief Online RLBWT construction.
 * @author Tomohiro I
 * @date 2017-10-12
 */
#ifndef INCLUDE_GUARD_OnlineRlbwt
#define INCLUDE_GUARD_OnlineRlbwt

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <fstream>

namespace itmmti
{
    using bwtintvl = std::pair<uint64_t, uint64_t>;

    /*!
     * @brief Online Run-length encoded Burrows–Wheeler transform (RLBWT).
     * @note
     *   ::OnlineRlbwt wraps ::DynRle to use it for representing dynamic RLE of BWT.
     *   In contrast to ::DynRle, OnlineRlbwt has an implicit end marker (em_) at emPos_.
     */
    template <class DynRle>
    class OnlineRlbwt
    {
    public:
        using CharT = typename DynRle::CharT;
        using BTreeNodeT = typename DynRle::BTreeNodeT;
        static constexpr uintptr_t NOTFOUND{BTreeNodeT::NOTFOUND};
        static constexpr uint8_t kB{DynRle::kB};
        static constexpr uint8_t kBtmBM{DynRle::kBtmBM};
        static constexpr uint8_t kBtmBS{DynRle::kBtmBS};

    private:
        DynRle drle_;
        uint64_t emPos_;       //!< Current position (0base) of end marker.
        CharT em_;             //!< End marker. It is used only when bwt[emPos_] is accessed (and does not matter if em_ appears in the input text).
        uint64_t num_em_;      // the number of em_
        uint64_t sap_s, sap_e; // cur sap interval [sap_s,sap_e]

    public:
        OnlineRlbwt(
            const size_t initNumBtms, //!< Initial size of DynRle to reserve.
            CharT em = 1              //!< End marker (default UINT64_MAX).
            ) : drle_(initNumBtms, 0),
                emPos_(0),
                em_(em)
        {
            num_em_ = 1;
            sap_s = 0;
            sap_e = 0;
        }

        /*!
         * @brief Get end marker.
         */
        CharT getEm() const noexcept
        {
            return em_;
        }

        /*!
         * @brief Get current position of end marker.
         */
        uint64_t getEndmarkerPos() const noexcept
        {
            return emPos_;
        }

        /*!
         * @brief Return current length including end marker.
         */
        uint64_t getLenWithEndmarker() const noexcept
        {
            return drle_.getSumOfWeight() + 1; // +1 for end marker, which is not in drle_.
        }

        /*!
         * @brief Extend RLBWT by appending one character.
         */
        void extend(
            const CharT ch //!< Character to append.
        )
        {
            uint64_t idxM = drle_.insertRun(emPos_, ch);
            emPos_ = drle_.rank(ch, idxM, emPos_, true);
        }

        void sptExtend(
            const CharT ch // !< Character to append.
        )
        {
            // std::cout << "(sap_s: " << sap_s << " - sap_e: " << sap_e << ")" << std::endl;
            // 插入当前元素
            if (sap_s == sap_e) // 代表没有后缀与目标后缀完全相同,只有一个地方可以插入
            {
                uint64_t tmp_sap_s = sap_s;
                drle_.insertRun(tmp_sap_s, ch);
            }
            else
            {
                // 计算有趣区间内是否有插入的符号
                uint64_t s_n = (sap_s == 0 ? 0 : drle_.rank(ch, sap_s - 1, false));
                uint64_t e_n = drle_.rank(ch, sap_e, false);
                // std::cout << "(s_n: " << s_n << " - e_n: " << e_n << ")" << std::endl;
                if (e_n - s_n > 0) // 有趣区间内已经存在插入的那个符号
                {
                    auto pos = drle_.select(ch, s_n + 1);
                    drle_.insertRun(pos, ch);
                }
                else // 有趣区间内不存在插入的那个符号
                {
                    // 插入原则是不从中打断别的run
                    drle_.optInsert(sap_s, sap_e, ch);
                }
            }
            // 计算下一个有趣区间
            if (ch == em_)
            {
                num_em_ += 1;
                sap_s = 0;
                sap_e = num_em_ - 1;
            }
            else
            {
                if (sap_s == sap_e)
                {
                    auto tmp = drle_.rank(ch, sap_s, true);
                    sap_s = tmp;
                    sap_e = tmp;
                }
                else
                {
                    sap_s = drle_.rank(ch, sap_s - 1, true) + 1;
                    sap_e = drle_.rank(ch, sap_e, true);
                }
            }
        }

        /*!
         * @brief Access to the current RLBWT by [] operator.
         */
        CharT operator[](
            uint64_t pos //!< in [0..OnlineRlbwt::getLenWithEndmarker()].
        ) const noexcept
        {
            assert(pos < getLenWithEndmarker());

            if (pos == emPos_)
            {
                return em_;
            }
            pos -= (pos > emPos_);
            uint64_t idxM = drle_.searchPosM(pos);
            return drle_.getCharFromIdxM(idxM);
        }

        /*!
         * @brief Return 'rank of ch at pos' + 'num of total occ of characters smaller than ch'.
         */
        uint64_t totalRank(
            const CharT ch,
            uint64_t pos //!< in [0..OnlineRlbwt::getLenWithEndmarker()].
        ) const noexcept
        {
            assert(pos < getLenWithEndmarker());

            pos -= (pos > emPos_);
            return drle_.rank(ch, pos, true);
        }

        /*!
         * @brief Compute bwt-interval for cW from bwt-interval for W
         * @note Intervals are [left, right) : right bound is excluded
         */
        bwtintvl lfMap(
            const bwtintvl intvl,
            const CharT ch) const noexcept
        {
            assert(intvl.first <= getLenWithEndmarker() && intvl.second <= getLenWithEndmarker());

            //// If "ch" is not in the alphabet or empty interval, return empty interval.
            const auto *retRootS = drle_.searchCharA(ch);
            if (retRootS->isDummy() || drle_.getCharFromNodeS(retRootS) != ch || intvl.first >= intvl.second)
            {
                return {0, 0};
            }

            uint64_t l = intvl.first - (intvl.first > emPos_);
            uint64_t r = intvl.second - (intvl.second > emPos_);
            const uint64_t idxM = drle_.searchPosM(l); // l is modified to relative pos.
            //// +1 because in F.select(0, ch) we are not taking into account the end-marker,
            //// which is in position 0 but not explicitly stored in F.
            return {
                drle_.rank(ch, idxM, l, true) - (drle_.getCharFromIdxM(idxM) == ch) + 1,
                drle_.rank(ch, r - 1, true) + 1};
        }

        /*!
         * @brief LF map.
         */
        uint64_t lfMap(uint64_t i)
        {
            assert(i < getLenWithEndmarker());

            i -= (i > emPos_);
            const uint64_t idxM = drle_.searchPosM(i);
            const auto ch = drle_.getCharFromIdxM(idxM);
            return drle_.rank(ch, idxM, i, true);
        }

        /*!
         * @brief Output original text to std::ofstream.
         */
        void invert(
            std::ofstream &ofs) const noexcept
        {
            uint64_t pos = 0;
            for (uint64_t i = 0; i < this->getLenWithEndmarker() - 1; ++i)
            {
                pos -= (pos > emPos_);
                const uint64_t idxM = drle_.searchPosM(pos);
                const auto ch = drle_.getCharFromIdxM(idxM);
                ofs.put(static_cast<signed char>(static_cast<unsigned char>(ch)));
                pos = drle_.rank(ch, idxM, pos, true);
            }
        }

        //////////////////////////////// statistics
        /*!
         * @brief Calculate total memory usage in bytes.
         */
        size_t calcMemBytes(
            bool includeThis = true) const noexcept
        {
            size_t size = sizeof(*this) * includeThis;
            size += drle_.calcMemBytes();
            return size;
        }

        /*!
         * @brief Print statistics
         */
        void printStatistics(
            std::ostream &os, //!< std::ostream (e.g., std::cout).
            const bool verbose) const noexcept
        {
            os << "OnlineRlbwt object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
            os << "Len with endmarker = " << getLenWithEndmarker() << std::endl;
            os << "emPos_ = " << emPos_ << ", em_ = " << em_ << std::endl;
            drle_.printStatistics(os, verbose);
            os << "OnlineRlbwt object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
        }

        void printDebugInfo(
            std::ostream &os //!< std::ostream (e.g., std::cout).
        ) const noexcept
        {
            os << "Len with endmarker = " << getLenWithEndmarker() << std::endl;
            os << "emPos_ = " << emPos_ << ", em_ = " << em_ << std::endl;
            drle_.printDebugInfo(os);
        }

        void printDetailInfo()
        {
            std::cout << "---------- bwt -----------" << std::endl;
            // std::cout << drle_.rank('p', 5, true) << std::endl;
            this->drle_.printDetailInfo();
            std::cout << std::endl;
        }

        void writeBWT(std::ofstream &ofs)
        {
            this->drle_.printString(ofs);
        }

        bool checkDecompress(
            std::ifstream &ifs) const noexcept
        {
            // {//debug
            //   std::cerr << __func__ << std::endl;
            // }

            uint64_t pos = 0;
            for (uint64_t i = 0; i < this->getLenWithEndmarker() - 1; ++i)
            {
                pos -= (pos > emPos_);
                // std::cout << "T[" << i << "] searchPos(" << pos << ") ";
                const uint64_t idxM = drle_.searchPosM(pos);
                const auto ch = static_cast<unsigned char>(drle_.getCharFromIdxM(idxM));
                // std::cout << "idxM = " << idxM << ", pos = " << pos << ", ch = " << (int)ch << "(" << ch << ")" << std::endl;
                char c; // Assume that the input character fits in char.
                ifs.get(c);
                unsigned char uc = static_cast<unsigned char>(c);
                if (uc != ch)
                {
                    std::cerr << "error: bad expansion at i = " << i << ", (" << ch << ") should be (" << uc << ")" << ", idxM = " << idxM << ", pos = " << pos << std::endl;
                    return false;
                }
                pos = drle_.rank(ch, idxM, pos, true);
            }
            return true;
        }
    };
};

#endif
