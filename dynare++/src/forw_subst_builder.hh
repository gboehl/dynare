/*
 * Copyright © 2006 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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

#ifndef FORW_SUBST_BUILDER_H
#define FORW_SUBST_BUILDER_H

#include <map>

#include "dynare_atoms.hh"

namespace ogdyn
{
  /* This struct encapsulates information about the process of forward
     substitutions. */
  struct ForwSubstInfo
  {
    int num_affected_equations{0};
    int num_subst_terms{0};
    int num_aux_variables{0};
    int num_new_terms{0};
  };

  class DynareModel;

  class ForwSubstBuilder
  {
    using Ttermauxmap = map<int, string>;
  protected:
    /* Reference to the model, to which we will add equations and change some
       equations. */
    DynareModel &model;
    /* A map mapping new auxiliary variables to the terms in the tree in the
       DynareModel. */
    Tsubstmap aux_map;
    /* Information about the substitutions. */
    ForwSubstInfo info;
  public:
    /* Do all the jobs needed. This scans all equations in the model, and for
       equations containing forward looking variables greater than 1 lead, it
       makes corresponding substitutions. Basically, it breaks each equation to
       its non-linear components and creates substitutions for these
       components, not for whole equation. This is because the expectation
       operator can go through the linear part of the function. This will save
       us many occurrences of other variables involved in the equation. */
    ForwSubstBuilder(DynareModel &m);
    ForwSubstBuilder(const ForwSubstBuilder &b) = delete;
    /* Copy constructor with a new instance of the model. */
    ForwSubstBuilder(const ForwSubstBuilder &b, DynareModel &m);
    /* Return the auxiliary variable mapping. */
    const Tsubstmap &
    get_aux_map() const
    {
      return aux_map;
    }
    /* Return the information. */
    const ForwSubstInfo &
    get_info() const
    {
      return info;
    }
  private:
    /* This method takes a nonlinear term t, and if it has leads of greater
       than 1, then it substitutes the term for the new variable (or string of
       variables). Note that the substitution is done by
       DynamicAtoms::assign_variable. This means that the substitution is made
       for all other ocurrences of t in the model. So there is no need of
       tracking already substituted terms. The other two parameters are just
       for identification of the new auxiliary variables. When called from the
       constructor, i is an equation number, j is an order of the non-linear
       term in the equation. */
    void substitute_for_term(int t, int i, int j);
    /* This is called just at the end of the job. It unassigns all nulary terms
       with a lead greater than 1. */
    void unassign_gt_1_leads();
    /* This unassigns all leads greater than 1 of the given name. */
    void unassign_gt_1_leads(const string &name);
  };
};

#endif
