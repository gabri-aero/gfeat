#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_

#include <iostream>
#include <mutex>

enum class Verbosity { Silent = 0, Info = 1 };

class Logger {
  private:
    Verbosity verbosity;
    std::mutex mutex_;
    Logger() { this->set_verbosity(Verbosity::Info); };

  public:
    static Logger &instance() {
        static Logger logger;
        return logger;
    }

    void set_verbosity(Verbosity verbosity) { this->verbosity = verbosity; }
    Verbosity get_verbosity() const { return this->verbosity; }

    template <typename T> Logger &operator<<(const T &value) {
        if (verbosity != Verbosity::Silent) {
            std::lock_guard<std::mutex> lock(mutex_);
            std::cout << value;
        }
        return *this;
    }

    // Overload for manipulators
    Logger &operator<<(std::ostream &(*manip)(std::ostream &)) {
        if (verbosity != Verbosity::Silent) {
            std::lock_guard<std::mutex> lock(mutex_);
            std::cout << manip;
        }
        return *this;
    }
};

#endif // _LOGGER_HPP_