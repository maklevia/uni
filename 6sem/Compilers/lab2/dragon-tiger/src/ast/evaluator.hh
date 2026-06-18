#ifndef EVALUATOR_HH
#define EVALUATOR_HH

#include "nodes.hh"

namespace ast {
namespace eval {

class Evaluator : public types::ConstASTIntVisitor {
public:
  virtual int32_t visit(const types::IntegerLiteral &) override;
  virtual int32_t visit(const types::StringLiteral &) override;
  virtual int32_t visit(const types::BinaryOperator &) override;
  virtual int32_t visit(const types::Sequence &) override;
  virtual int32_t visit(const types::Let &) override;
  virtual int32_t visit(const types::Identifier &) override;
  virtual int32_t visit(const types::IfThenElse &) override;
  virtual int32_t visit(const types::VarDecl &) override;
  virtual int32_t visit(const types::FunDecl &) override;
  virtual int32_t visit(const types::FunCall &) override;
  virtual int32_t visit(const types::WhileLoop &) override;
  virtual int32_t visit(const types::ForLoop &) override;
  virtual int32_t visit(const types::Break &) override;
  virtual int32_t visit(const types::Assign &) override;
};

} // namespace eval
} // namespace ast

#endif // EVALUATOR_HH
