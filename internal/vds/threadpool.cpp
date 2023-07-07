#include "threadpool.hpp"

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>


ThreadPool::ThreadPool(const std::size_t thread_count) noexcept (false)
    : threads{ thread_count }
{
    this->start();
}

// TODO safeguard against double starts
void ThreadPool::start() { 
    std::for_each(
        this->threads.begin(),
        this->threads.end(),
        [this](std::thread& thread) {
            thread = std::thread(&ThreadPool::thread_loop, this);
        }
    );
}

void ThreadPool::stop() noexcept (true) {
    {
        std::unique_lock< std::mutex > lock(this->mutex);
        this->should_terminate = true;
    }

    this->condition.notify_all();

    std::for_each(
        this->threads.begin(),
        this->threads.end(),
        [](std::thread& thread) { thread.join(); }
    );

    this->threads.clear();
}

bool ThreadPool::busy() noexcept (true) {
    std::unique_lock< std::mutex > lock(this->mutex);
    return !this->tasks.empty() or this->running_tasks;
}



std::future< void > ThreadPool::enqueue(
    std::function< void() > func
) noexcept (true) {
    auto promise = std::make_shared< std::promise< void > >();
    auto future  = promise->get_future();

    std::function< void() > task = [promise, func]() {
        try {
            func();
        } catch(...) {
            promise->set_exception( std::current_exception() );
        }
    };

    {
        std::unique_lock< std::mutex > lock(this->mutex);
        this->tasks.push( std::move(task) );
    }

    this->condition.notify_one();

    return future;
}

std::size_t ThreadPool::thread_count() const noexcept (true) {
    return this->threads.size();
}

void ThreadPool::thread_loop() {
    while (true) {
        std::function< void() > task;
        {
            std::unique_lock< std::mutex > lock(this->mutex);
            this->condition.wait(lock, [this] {
                return !this->tasks.empty() or this->should_terminate;
            });

            if (this->should_terminate)
                return;

            task = this->tasks.front();
            this->tasks.pop();
            ++this->running_tasks;

        }

        task();

        {
            std::unique_lock< std::mutex > lock(this->mutex);
            --this->running_tasks;
        }
    }
}

ThreadPool::~ThreadPool() {
    this->stop();
}
