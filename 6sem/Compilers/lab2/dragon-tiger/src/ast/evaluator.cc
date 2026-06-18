#include "evaluator.hh"
#include "../utils/errors.hh"

namespace ast {
namespace eval {

int32_t Evaluator::visit(const types::IntegerLiteral &literal) {
  return literal.value;
}

int32_t Evaluator::visit(const types::StringLiteral &) {
  utils::error("String literals are not supported in evaluator");
  return 0;
}

int32_t Evaluator::visit(const types::BinaryOperator &op) {
  int32_t left = op.get_left().accept(*this);
  int32_t right = op.get_right().accept(*this);
  switch (op.op) {
  case types::o_plus:
    return left + right;
  case types::o_minus:
    return left - right;
  case types::o_times:
    return left * right;
  case types::o_divide:
    if (right == 0)
      utils::error(op.loc, "division by zero");
    return left / right;
  case types::o_eq:
    return left == right;
  case types::o_neq:
    return left != right;
  case types::o_lt:
    return left < right;
  case types::o_le:
    return left <= right;
  case types::o_gt:
    return left > right;
  case types::o_ge:
    return left >= right;
  default:
    utils::error(op.loc, "unknown binary operator");
    return 0;
  }
}

int32_t Evaluator::visit(const types::Sequence &seq) {
  const auto &exprs = seq.get_exprs();
  if (exprs.empty()) {
    utils::error(seq.loc, "empty sequence");
    return 0;
  }
  int32_t result = 0;
  for (auto expr : exprs) {
    result = expr->accept(*this);
  }
  return result;
}

int32_t Evaluator::visit(const types::IfThenElse &ite) {
  int32_t cond = ite.get_condition().accept(*this);
  if (cond != 0) {
    return ite.get_then_part().accept(*this);
  } else {
    return ite.get_else_part().accept(*this);
  }
}

int32_t Evaluator::visit(const types::Let &) {
  utils::error("Let not supported");
  return 0;
}
int32_t Evaluator::visit(const types::Identifier &) {
  utils::error("Identifier not supported");
  return 0;
}
int32_t Evaluator::visit(const types::VarDecl &) {
  utils::error("VarDecl not supported");
  return 0;
}
int32_t Evaluator::visit(const types::FunDecl &) {
  utils::error("FunDecl not supported");
  return 0;
}
int32_t Evaluator::visit(const types::FunCall &) {
  utils::error("FunCall not supported");
  return 0;
}
int32_t Evaluator::visit(const types::WhileLoop &) {
  utils::error("WhileLoop not supported");
  return 0;
}
int32_t Evaluator::visit(const types::ForLoop &) {
  utils::error("ForLoop not supported");
  return 0;
}
int32_t Evaluator::visit(const types::Break &) {
  utils::error("Break not supported");
  return 0;
}
int32_t Evaluator::visit(const types::Assign &) {
  utils::error("Assign not supported");
  return 0;
}

} // namespace eval
} // namespace ast
