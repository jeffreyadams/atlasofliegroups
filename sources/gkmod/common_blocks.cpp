/*
  This is common_blocks.cpp

  Copyright (C) 2019,2020 Marc van Leeuwen
  Part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/


/*
  Motivation for this compilation unit.

  The computation of KLV polynomials is tightly bound to the notion of block,
  and indeed in takes place in a |kl::KL_table| value that is constructed from
  a block, and owned by a pointer member |kl_tab_ptr| of |blocks::Block_base|.

  In unitarity computations, blocks are computed from |Param| values, which
  internally correspond to the type |StandardRepr|; each one gives a single
  element of a block. This implies each time an entire block, or at least the
  part below the intitial parameter, is generated. To avoid excessive
  wastefulness, values derived from the KLV polynomials are associated in a
  permanent way with the parameters of the block in a |repr::Rep_table|, so that
  for a parameter whose block has already been subject to KLV computations, a
  second computation is avoided. However many parameters can be seen to have
  isomorphic blocks, namely when they share their |KGB| set and the integral
  subdatum (as determined by the infinitesimal character |gamma|, indeed already
  by its class modulo $X^*$); in this unit, blocks of type |common_block| are
  constructed that only use such information, and still allow computation of KLV
  polynomials (ordinary and twisted).
 */

namespace atlas {

} // |namespace atlas|
