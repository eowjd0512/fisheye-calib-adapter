#include "model/EUCM.hpp"

namespace FCA
{
namespace model
{

EUCM::EUCM(const std::string & model_name, const std::string & config_path)
: Base(model_name, config_path)
{
}

void EUCM::parse(){
    
}
void EUCM::initialize() {}
void EUCM::optimize() {}
Eigen::Vector2d EUCM::project() const {}
Eigen::Vector3d EUCM::unproject() const {}
Eigen::MatrixXd EUCM::calculate_jacobian() const {}
void EUCM::print() const {}
}  // namespace model
}  // namespace FCA
