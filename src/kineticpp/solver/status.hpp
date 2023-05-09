// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace kineticpp {
namespace solver {


enum struct ErrorCode : int {
    success = 0,
    iterations = -1,
    decomposition = -2,
};


static constexpr auto explain(ErrorCode code) {
    switch (code) {
    case ErrorCode::success:
        return "Success";
    case ErrorCode::iterations:
        return "Maximum iterations exceeded";
    case ErrorCode::decomposition:
        return "Matrix decomposition failed";
    default:
        return "Unknown error code";
    }
}


}  // namespace solver
}  // namespace kineticpp