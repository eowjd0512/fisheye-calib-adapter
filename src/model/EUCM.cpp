#include "model/EUCM.hpp"

namespace FCA
{
namespace model
{

EUCM::EUCM(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void EUCM::print() const {}
void EUCM::project() const {}
void EUCM::unproject() const {}
void EUCM::calculate_jacobian() const {}
void EUCM::get_jacobian() const {}
}  // namespace model
}  // namespace FCA
