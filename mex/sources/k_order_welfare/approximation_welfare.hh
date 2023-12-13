/*
 * Copyright © 2021-2023 Dynare Team
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

#ifndef APPROXIMATION_WELFARE_HH
#define APPROXIMATION_WELFARE_HH

#include "journal.hh"
#include "k_ord_objective.hh"

#include <memory>

class ApproximationWelfare
{
  KordwDynare& welfare;
  double discount_factor;
  std::unique_ptr<FGSContainer> rule_ders;
  std::unique_ptr<FGSContainer> rule_ders_s;
  std::unique_ptr<FGSContainer> cond_ders;
  std::unique_ptr<FoldDecisionRule> cond_fdr;
  IntSequence nvs;
  Journal& journal;

public:
  ApproximationWelfare(KordwDynare& w, double discount_factor, const FGSContainer& rule_ders,
                       const FGSContainer& rule_ders_s, Journal& j);

  [[nodiscard]] const KordwDynare&
  getWelfare() const
  {
    return welfare;
  }
  [[nodiscard]] const FGSContainer&
  get_rule_ders() const
  {
    return *rule_ders;
  }
  [[nodiscard]] const FGSContainer&
  get_rule_ders_s() const
  {
    return *rule_ders_s;
  }
  [[nodiscard]] const FGSContainer&
  get_cond_ders() const
  {
    return *cond_ders;
  }
  [[nodiscard]] const FoldDecisionRule& getFoldCondWel() const;

  void approxAtSteady();

protected:
  void saveRuleDerivs(const FGSContainer& W);
};

#endif
