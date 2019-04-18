// The contents of this file may be used under the terms of the Apache License,
// Version 2.0.
//
//    (See accompanying file LICENSE-Apache or copy at
//     http://www.apache.org/licenses/LICENSE-2.0)
//
// Alternatively, the contents of this file may be used under the terms of
// the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE-Boost or copy at
//     https://www.boost.org/LICENSE_1_0.txt)
//
// Unless required by applicable law or agreed to in writing, this software
// is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.

#ifndef RYU_CONFIG_H
#define RYU_CONFIG_H


#if defined(RYU_HEADER_ONLY)

#  if ! defined(__cplusplus)
#    error Header-only mode only available in C++
#  endif

#  define RYU_CPP_API

#  define RYU_INLINE inline
#  define RYU_PRIVATE_FUNC inline
#  define RYU_PUBLIC_FUNC inline

#  define RYU_EXTERN_C_BEGIN
#  define RYU_EXTERN_C_END

#else // defined(RYU_HEADER_ONLY)

#  if defined(__cplusplus)
#    define RYU_EXTERN_C_BEGIN extern "C" {
#    define RYU_EXTERN_C_END }
#  else
#    define RYU_EXTERN_C_BEGIN
#    define RYU_EXTERN_C_END
#  endif //defined(__cplusplus)

#  define RYU_INLINE static inline
#  define RYU_PRIVATE_FUNC static
#  define RYU_PUBLIC_FUNC extern

#endif // defined(RYU_HEADER_ONLY)

#if defined(RYU_CPP_API)
#  if ! defined(__cplusplus)
#    error C++ Header included from a non C++ source
#  endif
#  if ! defined(RYU_NAMESPACE_BEGIN)
#    define RYU_NAMESPACE_BEGIN namespace ryu {
#    define RYU_NAMESPACE_END }
#  elif ! defined(RYU_NAMESPACE_END)
#    error If you define RYU_NAMESPACE_BEGIN then you must define RYU_NAMESPACE_END
#  endif
#  define RYU_NAMESPACE_DETAIL_BEGIN namespace detail {
#  define RYU_NAMESPACE_DETAIL_END }
#  define RYU_USING_NAMESPACE_DETAIL using namespace detail;
#  define RYU_USING_NAMESPACE using namespace ryu;
#else // defined(RYU_CPP_API)
#  define RYU_NAMESPACE_BEGIN
#  define RYU_NAMESPACE_END
#  define RYU_NAMESPACE_DETAIL_BEGIN
#  define RYU_NAMESPACE_DETAIL_END
#  define RYU_USING_NAMESPACE_DETAIL
#  define RYU_USING_NAMESPACE
#endif // defined(RYU_CPP_API)

RYU_NAMESPACE_BEGIN
RYU_NAMESPACE_DETAIL_BEGIN
RYU_NAMESPACE_DETAIL_END
RYU_NAMESPACE_END

#endif // RYU_CONFIG_H
