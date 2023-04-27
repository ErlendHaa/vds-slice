#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>
#include <stdexcept>

#include "attribute.hpp"

Horizon::HorizontalIt Horizon::begin() const noexcept (true) {
    return HorizontalIt(this->m_ptr, this->vertical_window().size());
}

Horizon::HorizontalIt Horizon::end() const noexcept (true) {
    return HorizontalIt(
        this->m_ptr + this->size(),
        this->vertical_window().size()
    );
}
