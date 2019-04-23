// Copyright © 2006-2011, Ondra Kamenik

#include "forw_subst_builder.hh"

#include "dynare_model.hh"

using namespace ogdyn;

ForwSubstBuilder::ForwSubstBuilder(DynareModel &m)
  : model(m)
{
  info.num_new_terms -= model.getParser().getTree().get_num_op();

  // go through all equations
  int neq = model.eqs.nformulas();
  for (int i = 0; i < neq; i++)
    {
      int ft = model.eqs.formula(i);
      int mlead, mlag;
      model.termspan(ft, mlead, mlag);
      // if equation is too forward looking
      if (mlead > 1)
        {
          info.num_affected_equations++;
          // break it to non-linear terms
          unordered_set<int> nlt = model.get_nonlinear_subterms(ft);
          int j = 0; // indexes subterms
          // and make substitutions for all these non-linear subterms
          for (const auto &it : nlt)
            substitute_for_term(it, i, j++);
        }
    }
  // unassign all variables with lead greater than 1
  unassign_gt_1_leads();

  // forget the derivatives in the tree because some variables could
  // have been unassigned
  model.eqs.getTree().forget_derivative_maps();

  info.num_new_terms += model.getParser().getTree().get_num_op();
}

void
ForwSubstBuilder::substitute_for_term(int t, int i, int j)
{
  int mlead, mlag;
  model.termspan(t, mlead, mlag);
  if (mlead > 1)
    {
      info.num_subst_terms++;
      // Example for comments: let t = f(x(+4))
      // first make lagsubst be substitution setting f(x(+4)) to f(x(+1))
      // this is lag = -3 (1-mlead)
      map<int, int> lagsubst;
      unordered_set<int> nult = model.eqs.nulary_of_term(t);// make copy of nult!
      model.variable_shift_map(nult, 1-mlead, lagsubst);
      int lagt = model.eqs.add_substitution(t, lagsubst);

      // now maxlead of lagt is +1
      // add AUXLD_*_*_1 = f(x(+1)) to the model
      std::string name = "AUXLD_" + std::to_string(i) + '_' + std::to_string(j) + "_1";
      model.atoms.register_uniq_endo(name.c_str());
      info.num_aux_variables++;
      const char *ss = model.atoms.get_name_storage().query(name.c_str());
      int auxt = model.eqs.add_nulary(name.c_str());
      model.eqs.add_formula(model.eqs.add_binary(ogp::code_t::MINUS, auxt, lagt));
      aux_map.emplace(ss, lagt);
      // now add variables and equations
      // AUXLD_*_*_2 = AUXLD_*_*_1(+1) through
      // AUXLD_*_*_{mlead-1} = AUXLD_*_*_{mlead-2}(+1)
      for (int ll = 1; ll <= mlead-2; ll++)
        {
          // create AUXLD_*_*_{ll}(+1)
          name = "AUXLD_" + std::to_string(i) + '_' + std::to_string(j) + '_' + std::to_string(ll) + "(+1)";
          int lastauxt_lead = model.eqs.add_nulary(name.c_str());
          // create AUXLD_*_*{ll+1}
          name = "AUXLD_" + std::to_string(i) + '_' + std::to_string(j) + '_' + std::to_string(ll+1);
          model.atoms.register_uniq_endo(name.c_str());
          info.num_aux_variables++;
          ss = model.atoms.get_name_storage().query(name.c_str());
          auxt = model.eqs.add_nulary(name.c_str());
          // add AUXLD_*_*_{ll+1} = AUXLD_*_*_{ll}(+1)
          model.eqs.add_formula(model.eqs.add_binary(ogp::code_t::MINUS, auxt, lastauxt_lead));
          // add substitution to the map; todo: this
          // works well because in the context where
          // aux_map is used the timing doesn't matter,
          // however, it is misleading, needs to be
          // changed
          aux_map.emplace(ss, lagt);
        }

      // now we have to substitute AUXLD_*_*{mlead-1}(+1) for t
      name = "AUXLD_" + std::to_string(i) + '_' + std::to_string(j) + '_' + std::to_string(mlead-1);
      ss = model.atoms.get_name_storage().query(name.c_str());
      model.substitute_atom_for_term(ss, +1, t);
    }
}

void
ForwSubstBuilder::unassign_gt_1_leads(const char *name)
{
  const char *ss = model.atoms.get_name_storage().query(name);
  int mlead, mlag;
  model.atoms.varspan(name, mlead, mlag);
  for (int ll = 2; ll <= mlead; ll++)
    {
      int t = model.atoms.index(ss, ll);
      if (t != -1)
        model.atoms.unassign_variable(ss, ll, t);
    }
}

void
ForwSubstBuilder::unassign_gt_1_leads()
{
  const vector<const char *> &endovars = model.atoms.get_endovars();
  for (auto endovar : endovars)
    unassign_gt_1_leads(endovar);
  const vector<const char *> &exovars = model.atoms.get_exovars();
  for (auto exovar : exovars)
    unassign_gt_1_leads(exovar);
}

ForwSubstBuilder::ForwSubstBuilder(const ForwSubstBuilder &b, DynareModel &m)
  : model(m)
{
  for (auto it : b.aux_map)
    {
      const char *ss = m.atoms.get_name_storage().query(it.first);
      aux_map.emplace(ss, it.second);
    }
}
