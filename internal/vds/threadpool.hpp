#ifndef VDS_SLICE_TREADPOOL_HPP
#define VDS_SLICE_TREADPOOL_HPP

#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

struct Task {
    std::function< void() > func;
    std::promise< void > promise;
};

struct ThreadPool {
public:
    ThreadPool(const std::size_t nthreads) noexcept (false);

    ThreadPool(const ThreadPool& other) = delete;
    ThreadPool& operator=(const ThreadPool& other) = delete;

    void stop() noexcept (true);
    bool busy() noexcept (true);
    std::future< void > enqueue(std::function<void()> task) noexcept (true);
    std::size_t thread_count() const noexcept (true);

    ~ThreadPool();

private:
    void start();
    void thread_loop();

    bool should_terminate = false;
    std::size_t running_tasks;
    std::mutex mutex;
    std::condition_variable condition;

    std::vector< std::thread > threads;
    std::queue< std::function< void() > > tasks;
};

#endif // VDS_SLICE_TREADPOOL_HPP
