// Copyright © 2005-2011, Ondra Kamenik

#include "utils/cc/exception.hh"

#include "tree.hh"

#include <cmath>
#include <limits>

using namespace ogp;

/** Here we initialize OperationTree to contain only zero, one, nan
 * and two_over_pi terms. */
OperationTree::OperationTree()
{
  last_nulary = -1;
  // allocate space for the constants
  for (int i = 0; i < num_constants; i++)
    add_nulary();
}

int
OperationTree::add_nulary()
{
  int op = terms.size();
  terms.push_back({});
  _Tintset s;
  s.insert(op);
  nul_incidence.push_back(s);
  derivatives.push_back({});
  last_nulary = op;
  return op;
}

int
OperationTree::add_unary(code_t code, int op)
{
  if (op == zero
      && (code == code_t::UMINUS
          || code == code_t::SIN
          || code == code_t::TAN
          || code == code_t::SQRT
          || code == code_t::ERF))
    return zero;
  if ((op == zero && code == code_t::LOG) || op == nan)
    return nan;
  if (op == zero && (code == code_t::EXP
                     || code == code_t::COS
                     || code == code_t::ERFC))
    return one;

  Operation unary(code, op);
  auto i = opmap.find(unary);
  if (i == opmap.end())
    {
      int newop = terms.size();
      // add to the terms
      terms.push_back(unary);
      // copy incidence of the operand
      nul_incidence.push_back(nul_incidence[op]);
      // insert it to opmap
      opmap.emplace(unary, newop);
      // add empty map of derivatives
      _Tderivmap empty;
      derivatives.push_back(empty);
      return newop;
    }
  return i->second;
}

int
OperationTree::add_binary(code_t code, int op1, int op2)
{
  // quick exits for special values
  if (op1 == nan || op2 == nan)
    return nan;
  // for plus
  if (code == code_t::PLUS)
    {
      if (op1 == zero && op2 == zero)
        return zero;
      else if (op1 == zero)
        return op2;
      else if (op2 == zero)
        return op1;
    }
  // for minus
  if (code == code_t::MINUS)
    {
      if (op1 == zero && op2 == zero)
        return zero;
      else if (op1 == zero)
        return add_unary(code_t::UMINUS, op2);
      else if (op2 == zero)
        return op1;
    }
  // for times
  if (code == code_t::TIMES)
    {
      if (op1 == zero || op2 == zero)
        return zero;
      else if (op1 == one)
        return op2;
      else if (op2 == one)
        return op1;
    }
  // for divide
  if (code == code_t::DIVIDE)
    {
      if (op1 == op2)
        return one;
      else if (op1 == zero)
        return zero;
      else if (op2 == zero)
        return nan;
    }
  // for power
  if (code == code_t::POWER)
    {
      if (op1 == zero && op2 == zero)
        return nan;
      else if (op1 == zero)
        return zero;
      else if (op2 == zero)
        return one;
      else if (op1 == one)
        return one;
      else if (op2 == one)
        return op1;
    }

  // order operands of commutative operations
  if (code == code_t::TIMES || code == code_t::PLUS)
    if (op1 > op2)
      {
        int tmp = op1;
        op1 = op2;
        op2 = tmp;
      }

  // construct operation and check/add it
  Operation binary(code, op1, op2);
  auto i = opmap.find(binary);
  if (i == opmap.end())
    {
      int newop = terms.size();
      terms.push_back(binary);
      // sum both sets of incidenting nulary operations
      nul_incidence.push_back(nul_incidence[op1]);
      nul_incidence.back().insert(nul_incidence[op2].begin(), nul_incidence[op2].end());
      // add to opmap
      opmap.emplace(binary, newop);
      // add empty map of derivatives
      _Tderivmap empty;
      derivatives.push_back(empty);
      return newop;
    }
  return i->second;
}

int
OperationTree::add_derivative(int t, int v)
{
  if (t < 0 || t >= static_cast<int>(terms.size()))
    throw ogu::Exception(__FILE__, __LINE__,
                         "Wrong value for tree index in OperationTree::add_derivative");

  // quick returns for nulary terms or empty incidence
  if (terms[t].nary() == 0 && t != v)
    return zero;

  if (terms[t].nary() == 0 && t == v)
    return one;

  if (nul_incidence[t].end() == nul_incidence[t].find(v))
    return zero;

  // quick return if the derivative has been registered
  auto i = derivatives[t].find(v);
  if (i != derivatives[t].end())
    return i->second;

  int res = -1;
  switch (terms[t].getCode())
    {
    case code_t::UMINUS:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_unary(code_t::UMINUS, tmp);
        break;
      }
    case code_t::LOG:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_binary(code_t::DIVIDE, tmp, terms[t].getOp1());
        break;
      }
    case code_t::EXP:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_binary(code_t::TIMES, t, tmp);
        break;
      }
    case code_t::SIN:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_binary(code_t::TIMES, add_unary(code_t::COS, terms[t].getOp1()), tmp);
        break;
      }
    case code_t::COS:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_unary(code_t::UMINUS, add_binary(code_t::TIMES, add_unary(code_t::SIN, terms[t].getOp1()), tmp));
        break;
      }
    case code_t::TAN:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_unary(code_t::COS, terms[t].getOp1());
        res = add_binary(code_t::DIVIDE, tmp, add_binary(code_t::TIMES, tmp2, tmp2));
        break;
      }
    case code_t::SQRT:
      {
        int tmp = add_derivative(terms[t].getOp1(), v);
        res = add_binary(code_t::DIVIDE, tmp,
                         add_binary(code_t::PLUS, t, t));
        break;
      }
    case code_t::ERF:
      {
        int tmp = add_binary(code_t::TIMES, terms[t].getOp1(), terms[t].getOp1());
        tmp = add_unary(code_t::UMINUS, tmp);
        tmp = add_unary(code_t::EXP, tmp);
        int der = add_derivative(terms[t].getOp1(), v);
        tmp = add_binary(code_t::TIMES, tmp, der);
        res = add_binary(code_t::TIMES, two_over_pi, tmp);
        break;
      }
    case code_t::ERFC:
      {
        int tmp = add_binary(code_t::TIMES, terms[t].getOp1(), terms[t].getOp1());
        tmp = add_unary(code_t::UMINUS, tmp);
        tmp = add_unary(code_t::EXP, tmp);
        int der = add_derivative(terms[t].getOp1(), v);
        tmp = add_binary(code_t::TIMES, tmp, der);
        tmp = add_binary(code_t::TIMES, two_over_pi, tmp);
        res = add_unary(code_t::UMINUS, tmp);
        break;
      }
    case code_t::PLUS:
      {
        int tmp1 = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_derivative(terms[t].getOp2(), v);
        res = add_binary(code_t::PLUS, tmp1, tmp2);
        break;
      }
    case code_t::MINUS:
      {
        int tmp1 = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_derivative(terms[t].getOp2(), v);
        res = add_binary(code_t::MINUS, tmp1, tmp2);
        break;
      }
    case code_t::TIMES:
      {
        int tmp1 = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_derivative(terms[t].getOp2(), v);
        int res1 = add_binary(code_t::TIMES, terms[t].getOp1(), tmp2);
        int     res2 = add_binary(code_t::TIMES, tmp1, terms[t].getOp2());
        res = add_binary(code_t::PLUS, res1, res2);
        break;
      }
    case code_t::DIVIDE:
      {
        int tmp1 = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_derivative(terms[t].getOp2(), v);
        if (tmp2 == zero)
          res = add_binary(code_t::DIVIDE, tmp1, terms[t].getOp2());
        else
          {
            int nom = add_binary(code_t::MINUS,
                                 add_binary(code_t::TIMES, tmp1, terms[t].getOp2()),
                                 add_binary(code_t::TIMES, tmp2, terms[t].getOp1()));
            int den = add_binary(code_t::TIMES, terms[t].getOp2(), terms[t].getOp2());
            res = add_binary(code_t::DIVIDE, nom, den);
          }
        break;
      }
    case code_t::POWER:
      {
        int tmp1 = add_derivative(terms[t].getOp1(), v);
        int tmp2 = add_derivative(terms[t].getOp2(), v);
        int s1 = add_binary(code_t::TIMES, tmp2,
                            add_binary(code_t::TIMES, t,
                                       add_unary(code_t::LOG, terms[t].getOp1())));
        int s2 = add_binary(code_t::TIMES, tmp1,
                            add_binary(code_t::TIMES, terms[t].getOp2(),
                                       add_binary(code_t::POWER, terms[t].getOp1(),
                                                  add_binary(code_t::MINUS, terms[t].getOp2(), one))));
        res = add_binary(code_t::PLUS, s1, s2);
        break;
      }
    case code_t::NONE:
      break;
    }

  if (res == -1)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Unknown operation code.");

  register_derivative(t, v, res);

  return res;
}

int
OperationTree::add_substitution(int t, const map<int, int> &subst)
{
  return add_substitution(t, subst, *this);
}

int
OperationTree::add_substitution(int t, const map<int, int> &subst,
                                const OperationTree &otree)
{
  // return substitution of t if it is in the map
  auto it = subst.find(t);
  if (subst.end() != it)
    return it->second;

  int nary = otree.terms[t].nary();
  if (nary == 2)
    {
      // return the binary operation of the substituted terms
      int t1 = add_substitution(otree.terms[t].getOp1(), subst, otree);
      int t2 = add_substitution(otree.terms[t].getOp2(), subst, otree);
      return add_binary(otree.terms[t].getCode(), t1, t2);
    }
  else if (nary == 1)
    {
      // return the unary operation of the substituted term
      int t1 = add_substitution(otree.terms[t].getOp1(), subst, otree);
      return add_unary(otree.terms[t].getCode(), t1);
    }
  else
    {
      // if t is not the first num_constants, and otree is not this
      // tree, then raise and exception. Otherwise return t, since
      // it is either a special term (having the same semantics in
      // both trees), or the trees are the same, hence t has the
      // same semantics
      if (t < num_constants || this == &otree)
        return t;
      else
        {
          throw ogu::Exception(__FILE__, __LINE__,
                               "Incomplete substitution map in OperationTree::add_substitution");
          return -1;
        }
    }
}

void
OperationTree::nularify(int t)
{
  // remove the original operation from opmap
  auto it = opmap.find(terms[t]);
  if (it != opmap.end())
    opmap.erase(it);
  // turn the operation to nulary
  Operation nulary_op;
  terms[t] = nulary_op;
  // update last nulary
  if (last_nulary < t)
    last_nulary = t;
  // update nul_incidence information for all terms including t
  update_nul_incidence_after_nularify(t);
}

void
OperationTree::register_derivative(int t, int v, int tder)
{
  // todo: might check that the insert inserts a new pair
  derivatives[t].emplace(v, tder);
}

unordered_set<int>
OperationTree::select_terms(int t, const opselector &sel) const
{
  unordered_set<int> subterms;
  select_terms(t, sel, subterms);
  return subterms;
}

void
OperationTree::select_terms(int t, const opselector &sel, unordered_set<int> &subterms) const
{
  const Operation &op = terms[t];

  if (sel(t))
    subterms.insert(t);
  else
    if (op.nary() == 2)
      {
        select_terms(op.getOp1(), sel, subterms);
        select_terms(op.getOp2(), sel, subterms);
      }
    else if (op.nary() == 1)
      select_terms(op.getOp1(), sel, subterms);
}

unordered_set<int>
OperationTree::select_terms_inv(int t, const opselector &sel) const
{
  unordered_set<int> subterms;
  select_terms_inv(t, sel, subterms);
  return subterms;
}

bool
OperationTree::select_terms_inv(int t, const opselector &sel, unordered_set<int> &subterms) const
{
  const Operation &op = terms[t];

  if (op.nary() == 2)
    {
      bool a1 = select_terms_inv(op.getOp1(), sel, subterms);
      bool a2 = select_terms_inv(op.getOp2(), sel, subterms);
      if (a1 && a2 && sel(t))
        {
          subterms.insert(t);
          return true;
        }
    }
  else if (op.nary() == 1)
    {
      bool a1 = select_terms_inv(op.getOp1(), sel, subterms);
      if (a1 && sel(t))
        {
          subterms.insert(t);
          return true;
        }
    }
  else
    {
      if (sel(t))
        {
          subterms.insert(t);
          return true;
        }
    }

  return false;
}

void
OperationTree::forget_derivative_maps()
{
  for (auto & derivative : derivatives)
    derivative.clear();
}

void
OperationTree::print_operation_tree(int t, std::ostream &os, OperationFormatter &f) const
{
  f.format(terms[t], t, os);
}

void
OperationTree::print_operation(int t) const
{
  DefaultOperationFormatter dof(*this);
  print_operation_tree(t, std::cout, dof);
}

void
OperationTree::update_nul_incidence_after_nularify(int t)
{
  unordered_set<int> updated;
  for (int tnode = num_constants; tnode < static_cast<int>(terms.size()); tnode++)
    {
      const Operation &op = terms[tnode];
      if (op.nary() == 2)
        {
          int op1 = op.getOp1();
          int op2 = op.getOp2();
          if (op1 >= tnode || op2 >= tnode)
            throw ogu::Exception(__FILE__, __LINE__,
                                 "Tree disorder asserted");
          bool updated1 = (updated.end() != updated.find(op1));
          bool updated2 = (updated.end() != updated.find(op2));
          if (updated1 || updated2)
            {
              nul_incidence[tnode] = nul_incidence[op1];
              nul_incidence[tnode].insert(nul_incidence[op2].begin(), nul_incidence[op2].end());
              updated.insert(tnode);
            }
        }
      else if (op.nary() == 1)
        {
          int op1 = op.getOp1();
          if (op1 >= tnode)
            throw ogu::Exception(__FILE__, __LINE__,
                                 "Tree disorder asserted");
          bool updated1 = (updated.end() != updated.find(op1));
          if (updated1)
            {
              nul_incidence[tnode] = nul_incidence[op1];
              updated.insert(tnode);
            }
        }
      else if (op.nary() == 0)
        {
          if (tnode == t)
            {
              nul_incidence[tnode].clear();
              nul_incidence[tnode].insert(tnode);
              updated.insert(tnode);
            }
        }
    }
}

EvalTree::EvalTree(const OperationTree &ot, int last)
  : otree(ot),
    values(std::make_unique<double[]>((last == -1) ? ot.terms.size() : last+1)),
    flags(std::make_unique<bool[]>((last == -1) ? ot.terms.size() : last+1)),
    last_operation((last == -1) ? ot.terms.size()-1 : last)
{
  if (last_operation < OperationTree::num_constants-1
      || last_operation > static_cast<int>(ot.terms.size())-1)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Wrong last in EvalTree constructor.");

  values[0] = 0.0;
  flags[0] = true;
  values[1] = 1.0;
  flags[1] = true;
  values[2] = std::numeric_limits<double>::quiet_NaN();
  flags[2] = true;
  values[3] = 2.0/sqrt(M_PI);
  flags[3] = true;
  // this sets from num_constants on
  reset_all();
}

void
EvalTree::reset_all()
{
  for (int i = OperationTree::num_constants; i <= last_operation; i++)
    flags[i] = false;
}

void
EvalTree::set_nulary(int t, double val)
{
  if (t < 0 || t > last_operation)
    throw ogu::Exception(__FILE__, __LINE__,
                         "The tree index out of bounds in EvalTree::set_nulary");
  if (t < OperationTree::num_constants || otree.terms[t].nary() != 0)
    throw ogu::Exception(__FILE__, __LINE__,
                         "The term is not nulary assignable in EvalTree::set_nulary");

  values[t] = val;
  flags[t] = true;
}

double
EvalTree::eval(int t)
{
  if (t < 0 || t > last_operation)
    throw ogu::Exception(__FILE__, __LINE__,
                         "The tree index out of bounds in EvalTree::eval");
  if (otree.terms[t].nary() == 0 && flags[t] == false)
    throw ogu::Exception(__FILE__, __LINE__,
                         "Nulary term has not been assigned a value in EvalTree::eval");

  if (!flags[t])
    {
      const Operation &op = otree.terms[t];
      if (op.nary() == 1)
        {
          double r1 = eval(op.getOp1());
          double res;
          if (op.getCode() == code_t::UMINUS)
            res = -r1;
          else if (op.getCode() == code_t::LOG)
            res = log(r1);
          else if (op.getCode() == code_t::EXP)
            res = exp(r1);
          else if (op.getCode() == code_t::SIN)
            res = sin(r1);
          else if (op.getCode() == code_t::COS)
            res = cos(r1);
          else if (op.getCode() == code_t::TAN)
            res = tan(r1);
          else if (op.getCode() == code_t::SQRT)
            res = sqrt(r1);
          else if (op.getCode() == code_t::ERF)
            res = erf(r1);
          else if (op.getCode() == code_t::ERFC)
            res = erfc(r1);
          else
            {
              throw ogu::Exception(__FILE__, __LINE__,
                                   "Unknown unary operation code in EvalTree::eval");
              res = 0.0;
            }
          values[t] = res;
          flags[t] = true;
        }
      else if (op.nary() == 2)
        {
          double res;
          if (op.getCode() == code_t::PLUS)
            {
              double r1 = eval(op.getOp1());
              double r2 = eval(op.getOp2());
              res = r1 + r2;
            }
          else if (op.getCode() == code_t::MINUS)
            {
              double r1 = eval(op.getOp1());
              double r2 = eval(op.getOp2());
              res = r1 - r2;
            }
          else if (op.getCode() == code_t::TIMES)
            {
              // pickup less complex formula first
              unsigned int nul1 = otree.nulary_of_term(op.getOp1()).size();
              unsigned int nul2 = otree.nulary_of_term(op.getOp2()).size();
              if (nul1 < nul2)
                {
                  double r1 = eval(op.getOp1());
                  if (r1 == 0.0)
                    res = 0.0;
                  else
                    {
                      double r2 = eval(op.getOp2());
                      res = r1 * r2;
                    }
                }
              else
                {
                  double r2 = eval(op.getOp2());
                  if (r2 == 0)
                    res = 0.0;
                  else
                    {
                      double r1 = eval(op.getOp1());
                      res = r1*r2;
                    }
                }
            }
          else if (op.getCode() == code_t::DIVIDE)
            {
              double r1 = eval(op.getOp1());
              if (r1 == 0)
                res = 0.0;
              else
                {
                  double r2 = eval(op.getOp2());
                  res = r1 / r2;
                }
            }
          else if (op.getCode() == code_t::POWER)
            {
              // suppose that more complex is the first op in average
              double r2 = eval(op.getOp2());
              if (r2 == 0.0)
                res = 1.0;
              else
                {
                  double r1 = eval(op.getOp1());
                  res = pow(r1, r2);
                }
            }
          else
            {
              throw ogu::Exception(__FILE__, __LINE__,
                                   "Unknown binary operation code in EvalTree::eval");
              res = 0.0;
            }
          values[t] = res;
          flags[t] = true;
        }
      return values[t];
    }

  // if (! std::isfinite(values[t]))
  //	printf("Tree value t=%d is not finite = %f\n", t, values[t]);

  return values[t];
}

void
EvalTree::print() const
{
  printf("last_op=%d\n", last_operation);
  printf("         0     1     2     3     4     5     6     7     8     9\n");
  printf("----------------------------------------------------------------\n");
  for (int i = 0; i <= (last_operation+1)/10; i++)
    {
      printf("%-3d|", i);
      int j = 0;
      while (j < 10 && 10*i+j < last_operation+1)
        {
          int k = 10*i+j;
          if (flags[k])
            printf(" %5.1g", values[k]);
          else
            printf(" -----");
          j++;
        }
      printf("\n");
    }
}

void
DefaultOperationFormatter::format(const Operation &op, int t, std::ostream &os)
{
  // add to the stop_set
  if (stop_set.end() == stop_set.find(t))
    stop_set.insert(t);
  else
    return;

  // call recursively non-nulary terms of the operation
  if (op.nary() == 2)
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.terms[t1];
      int t2 = op.getOp2();
      const Operation &op2 = otree.terms[t2];
      if (op1.nary() > 0)
        format(op1, t1, os);
      if (op2.nary() > 0)
        format(op2, t2, os);
    }
  if (op.nary() == 1)
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.terms[t1];
      if (op1.nary() > 0)
        format(op1, t1, os);
    }

  // print 'term ='
  format_term(t, os);
  os << " = ";
  if (op.nary() == 0)
    format_nulary(t, os);
  else if (op.nary() == 1)
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.terms[t1];
      std::string opname = "unknown";
      switch (op.getCode())
        {
        case code_t::UMINUS:
          opname = "-";
          break;
        case code_t::LOG:
          opname = "log";
          break;
        case code_t::EXP:
          opname = "exp";
          break;
        case code_t::SIN:
          opname = "sin";
          break;
        case code_t::COS:
          opname = "cos";
          break;
        case code_t::TAN:
          opname = "tan";
          break;
        case code_t::SQRT:
          opname = "sqrt";
          break;
        case code_t::ERF:
          opname = "erf";
          break;
        case code_t::ERFC:
          opname = "erfc";
          break;
        default:
          break;
        }
      os << opname << '(';
      if (op1.nary() == 0)
        format_nulary(t1, os);
      else
        format_term(t1, os);
      os << ")";
    }
  else
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.terms[t1];
      int t2 = op.getOp2();
      const Operation &op2 = otree.terms[t2];
      std::string opname = "unknown";
      switch (op.getCode())
        {
        case code_t::PLUS:
          opname = "+";
          break;
        case code_t::MINUS:
          opname = "-";
          break;
        case code_t::TIMES:
          opname = "*";
          break;
        case code_t::DIVIDE:
          opname = "/";
          break;
        case code_t::POWER:
          opname = "^";
          break;
        default:
          break;
        }
      if (op1.nary() == 0)
        format_nulary(t1, os);
      else
        format_term(t1, os);
      os << ' ' << opname << ' ';
      if (op2.nary() == 0)
        format_nulary(t2, os);
      else
        format_term(t2, os);
    }

  print_delim(os);
}

void
DefaultOperationFormatter::format_term(int t, std::ostream &os) const
{
  os << '$' << t;
}

void
DefaultOperationFormatter::format_nulary(int t, std::ostream &os) const
{
  if (t == OperationTree::zero)
    os << '0';
  else if (t == OperationTree::one)
    os << '1';
  else if (t == OperationTree::nan)
    os << "NaN";
  else
    os << '$' << t;
}

void
DefaultOperationFormatter::print_delim(std::ostream &os) const
{
  os << ";\n";
}

std::string
OperationStringConvertor::convert(const Operation &op, int t) const
{
  if (op.nary() == 0)
    {
      if (t < OperationTree::num_constants)
        if (t == OperationTree::zero)
          return std::string("0");
        else if (t == OperationTree::one)
          return std::string("1");
        else if (t == OperationTree::nan)
          return std::string("NaN");
        else if (t == OperationTree::two_over_pi)
          {
            char buf[100];
            sprintf(buf, "%20.16g", 2.0/std::sqrt(M_PI));
            return std::string(buf);
          }
        else
          {
            return std::string("error!error");
          }
      else
        return nulsc.convert(t);
    }
  else if (op.nary() == 1)
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.operation(t1);
      const char *opname = "unknown";
      switch (op.getCode())
        {
        case code_t::UMINUS:
          opname = "-";
          break;
        case code_t::LOG:
          opname = "log";
          break;
        case code_t::EXP:
          opname = "exp";
          break;
        case code_t::SIN:
          opname = "sin";
          break;
        case code_t::COS:
          opname = "cos";
          break;
        case code_t::TAN:
          opname = "tan";
          break;
        case code_t::SQRT:
          opname = "sqrt";
          break;
        case code_t::ERF:
          opname = "erf";
          break;
        case code_t::ERFC:
          opname = "erfc";
          break;
        default:
          break;
        }
      std::string s1 = convert(op1, t1);
      return std::string(opname) + "(" + s1 + ")";
    }
  else
    {
      int t1 = op.getOp1();
      const Operation &op1 = otree.operation(t1);
      int t2 = op.getOp2();
      const Operation &op2 = otree.operation(t2);
      const char *opname = "unknown";
      switch (op.getCode())
        {
        case code_t::PLUS:
          opname = "+";
          break;
        case code_t::MINUS:
          opname = "-";
          break;
        case code_t::TIMES:
          opname = "*";
          break;
        case code_t::DIVIDE:
          opname = "/";
          break;
        case code_t::POWER:
          opname = "^";
          break;
        default:
          break;
        }
      // decide about parenthesis
      bool op1_par = true;
      bool op2_par = true;
      if (op.getCode() == code_t::PLUS)
        {
          op1_par = false;
          op2_par = false;
        }
      else if (op.getCode() == code_t::MINUS)
        {
          op1_par = false;
          if (op2.getCode() != code_t::MINUS && op2.getCode() != code_t::PLUS)
            op2_par = false;
        }
      else
        {
          if (op1.nary() < 2)
            op1_par = false;
          if (op2.nary() < 2)
            op2_par = false;
        }

      std::string res;
      if (op1_par)
        res += "(";
      res += convert(op1, t1);
      if (op1_par)
        res += ")";
      res += " ";
      res += opname;
      res += " ";
      if (op2_par)
        res += "(";
      res += convert(op2, t2);
      if (op2_par)
        res += ")";

      return res;
    }
}
