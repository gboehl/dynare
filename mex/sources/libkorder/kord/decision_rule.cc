/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "kord_exception.hh"
#include "decision_rule.hh"
#include "dynamic_model.hh"

#include "SymSchurDecomp.hh"

#include <memory>

// FoldDecisionRule conversion from UnfoldDecisionRule
FoldDecisionRule::FoldDecisionRule(const UnfoldDecisionRule &udr)
  : DecisionRuleImpl<Storage::fold>(ctraits<Storage::fold>::Tpol(udr.nrows(), udr.nvars()),
                                    udr.ypart, udr.nu, udr.ysteady)
{
  for (const auto &it : udr)
    insert(std::make_unique<ctraits<Storage::fold>::Ttensym>(*(it.second)));
}

// UnfoldDecisionRule conversion from FoldDecisionRule
UnfoldDecisionRule::UnfoldDecisionRule(const FoldDecisionRule &fdr)
  : DecisionRuleImpl<Storage::unfold>(ctraits<Storage::unfold>::Tpol(fdr.nrows(), fdr.nvars()),
                                      fdr.ypart, fdr.nu, fdr.ysteady)
{
  for (const auto &it : fdr)
    insert(std::make_unique<ctraits<Storage::unfold>::Ttensym>(*(it.second)));
}
