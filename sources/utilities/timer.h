/*
  This is timer.h

  Copyright (C) 2023 Marc van Leeuwen
  part of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

/*
  The purpose of this module it to provide a simple class to record the time
  since the class instance was constructed. Since it must have a long lifetime
  in order to be useful, instances are most likely to be |static| variables.
*/

#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace atlas {
namespace time {

class Timer
{
  std::chrono::time_point<std::chrono::steady_clock> start;
public:
  Timer() : start(std::chrono::steady_clock::now()) {}
  long long int elapsed_ms() const // number of milliseconds since construction
  {
    using namespace std::literals::chrono_literals;
    return (std::chrono::steady_clock::now()-start)/1ms;
  }
}; // class |Timer|

} // |namespace time|
} // |namespace atlas|

#endif
