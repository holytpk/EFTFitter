// -*- C++ -*-
// author: afiq anuar
// short: common utility methods on data handling

#ifndef TEMPLATEUTIL_H
#define TEMPLATEUTIL_H

// headers are not necessarily needed HERE
// but convenient to have one central include that allows access to very common techniques

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <memory>
#include <utility>

#include <vector>
#include <map>
#include <array>
#include <tuple>

#include <algorithm>
#include <numeric>

#include <type_traits>
#include <functional>
#include <stdexcept>

#include <cmath>

/// poor man's std::hypot (tested with g++/clang++, fine for double: 1e-158 < i < 1e10)
template <typename Number>
Number quad_sum(const Number num1, const Number num2)
{
  return std::sqrt((num1 * num1) + (num2 * num2));
}



/// as it says on the tin
template <typename Number>
Number deltaPhi(const Number phi1, const Number phi2) 
{
  return std::acos( std::cos(phi1 - phi2) );
}



/// also as it says on the tin
template <typename Number>
Number deltaR(const Number eta1, const Number phi1, const Number eta2, const Number phi2) 
{
  return quad_sum(eta1 - eta2, deltaPhi(phi1, phi2));
}



/// string replacement
inline bool replace(std::string &str, const std::string &from, const std::string &to) {
  std::size_t start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}



/// count occurences of substrings
/// credit https://stackoverflow.com/questions/22489073/counting-the-number-of-occurrences-of-a-string-within-a-string
inline int count_substring(const std::string &str, const std::string &sub) {
  int nSub = 0;
  std::string::size_type iSub = 0;
  while ((iSub = str.find(sub, iSub)) != std::string::npos) {
    ++nSub;
    iSub += sub.length();
  }

  return nSub;
}



/// number to string; to_string tend to give more precision than needed
template <typename Number> 
std::string toStr(Number inNo, const int prec = -1, const bool fixed = false) 
{
  std::ostringstream outStr;
  if (fixed)
    outStr << std::fixed;
  if (prec > 0)
    outStr << std::setprecision(prec);

  outStr << inNo; 
  return outStr.str(); 
}



/// when filling vectors with intervals are needed; ensure last element in vec equals from
/// step serves also as epsilon (should be ok for intended usage of this function)
template <typename Number, typename std::enable_if_t<std::is_floating_point<Number>::value>* = nullptr> 
void fillInterval(std::vector<Number> &vec, Number from, Number to, Number step)
{
  // do nothing in various cases the function won't make sense
  const Number eps = std::abs(step / 50.);
  if (step == 0. or (step > 0. and to - from < step) or (step < 0. and to - from > step)
      or vec.empty() or std::abs(vec.back() - from) >= eps) return;

  while (std::abs(to - vec.back()) >= eps)
    vec.emplace_back(vec.back() + step);
}



/// also also for integral types where there's none of the epsilon shenanigans
template <typename Number, typename std::enable_if_t<std::is_integral<Number>::value>* = nullptr> 
void fillInterval(std::vector<Number> &vec, Number from, Number to, Number step)
{
  // do nothing in various cases the function won't make sense
  if (step == 0 or (step > 0 and to - from < step) or (step < 0 and to - from > step)
      or vec.empty() or vec.back() != from) return;

  while (std::abs(to - vec.back()) >= std::abs(step))
    vec.emplace_back(vec.back() + step);
}



/// method that makes such an interval vector relying on the above
template <typename Number> 
std::vector<Number> makeInterval(Number from, Number to, Number step)
{
  std::vector<Number> v_interval = {from};
  fillInterval(v_interval, from, to, step);
  return (v_interval.size() > 1) ? v_interval : std::vector<Number>{};
}



/// unary predicate that is always true; to be used in the specializations (where all elements are to be extracted)
/// usage of sizeof (which is never false) rather than literal true is to dodge -Werror=unused-parameter
const auto always_true = [] (const auto &obj) { return sizeof(obj); };



/// functions for returning pair elements
const auto return_first = [] (const auto &p) { return p.first; };
const auto return_second = [] (const auto &p) { return p.second; };



/// extract stuff from data contained in a map satisfying a predicate (or transformation thereof); output dumped in a vector
/// predicate check is done on map elements, before the transformation product (but one could nest them of course)
template <typename Key, typename Value, typename Function, typename Predicate> 
auto extractMapIf(const std::map<Key, Value> &map, Function func, Predicate pred)
{
  std::vector<decltype( func(*(std::cbegin(map))) )> v_out;
  for (const auto& p : map) {
    if (pred(p))
      v_out.emplace_back( func(p) );
  }

  return v_out;
}



/// specialization of the above where predicate is always true
template <typename Key, typename Value, typename Function> 
auto extractMap(const std::map<Key, Value> &map, Function func)
{
  return extractMapIf(map, func, always_true);
}



/// specialization of the above for extracting all keys
template <typename Key, typename Value> 
std::vector<Key> extractKey(const std::map<Key, Value> &map)
{
  return extractMapIf(map, return_first, always_true);
}



/// same as above, with a predicate
template <typename Key, typename Value, typename Predicate> 
std::vector<Key> extractKeyIf(const std::map<Key, Value> &map, Predicate pred)
{
  return extractMapIf(map, return_first, pred);
}



/// and now for value variant
template <typename Key, typename Value> 
std::vector<Value> extractValue(const std::map<Key, Value> &map)
{
  return extractMapIf(map, return_second, always_true);
}



/// and its predicate version
template <typename Key, typename Value, typename Predicate> 
std::vector<Value> extractValueIf(const std::map<Key, Value> &map, Predicate pred)
{
  return extractMapIf(map, return_second, pred);
}



/// this overload is for printAll to correctly print pair-types
template <typename One, typename Two> 
std::ostream& operator<<(std::ostream& out, const std::pair<One, Two> &pair)
{
  out << "[ " << pair.first << " :: " << pair.second << " ]";
  return out;
}



/// check if a given type is printable; don't really know how to use though - if (...) in printAll and/or << didn't work :(
/// credits Tony Delroy @ https://stackoverflow.com/questions/22758291/how-can-i-detect-if-a-type-can-be-streamed-to-an-stdostream
template <typename T> class is_streamable {
  // must be template to get SFINAE fall-through...
  template <typename U> static auto test(const U* u) -> decltype(std::cout << *u);
  static auto test(...) -> std::false_type;

  public:
    enum { value = !std::is_same<decltype(test( (T*) 0 )), std::false_type>::value };
};



/// master-cout of elements in container https://eli.thegreenplace.net/2014/variadic-templates-in-c/
/// it requires << to be defined for type Value however, so no printing of map of vectors etc yet :(
/// note to Afiq: don't, DON'T try to work on this; it's hard
/// wont compile if used for containers with not-streamable value types
template <template <typename, typename...> class Container, typename Value, typename... Args,
          typename std::enable_if_t<is_streamable< typename Container<Value, Args...>::value_type >::value>* = nullptr>
void printAll(const Container<Value, Args...> &con)
{
  for (const auto &ele : con)
    std::cout << ele << std::endl;
  std::cout << std::endl;
}

#endif
