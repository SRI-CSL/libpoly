#include <polyxx/context.h>

namespace poly {

  Context::Context() {
    mVariableDB = deleting_unique_ptr<lp_variable_db_t>(
        lp_variable_db_new(),
        [](lp_variable_db_t* ptr) { lp_variable_db_detach(ptr); });
    mVariableOrder = deleting_unique_ptr<lp_variable_order_t>(
        lp_variable_order_new(),
        [](lp_variable_order_t* ptr) { lp_variable_order_detach(ptr); });
    mPolynomialContext = deleting_unique_ptr<lp_polynomial_context_t>(
        lp_polynomial_context_new(lp_Z, mVariableDB.get(),
                                  mVariableOrder.get()),
        [](lp_polynomial_context_t* ptr) {
          lp_polynomial_context_detach(ptr);
        });
  }

  Context::Context(lp_polynomial_context_t* ctx) {
    mVariableDB = deleting_unique_ptr<lp_variable_db_t>(
      ctx->var_db,
      [](lp_variable_db_t *ptr) { lp_variable_db_detach(ptr); });
    mVariableOrder = deleting_unique_ptr<lp_variable_order_t>(
      ctx->var_order,
      [](lp_variable_order_t *ptr) { lp_variable_order_detach(ptr); });
    mPolynomialContext = deleting_unique_ptr<lp_polynomial_context_t>(
      ctx,
      [](lp_polynomial_context_t *ptr) { lp_polynomial_context_detach(ptr); });

    lp_variable_db_attach(ctx->var_db);
    lp_variable_order_attach(ctx->var_order);
    lp_polynomial_context_attach(ctx);
  }

  Context::Context(const Context& other) : Context(other.get_polynomial_context()) {}

  lp_variable_db_t* Context::get_variable_db() const {
    return const_cast<lp_variable_db_t*>(mVariableDB.get());
  }

  lp_variable_order_t* Context::get_variable_order() const {
    return const_cast<lp_variable_order_t*>(mVariableOrder.get());
  }

  lp_polynomial_context_t* Context::get_polynomial_context() const {
    return const_cast<lp_polynomial_context_t*>(mPolynomialContext.get());
  }

  const Context& Context::get_context() {
    static Context c;
    return c;
  }

}  // namespace poly