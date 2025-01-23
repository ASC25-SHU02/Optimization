[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11046885.svg)](https://doi.org/10.5281/zenodo.11046885)

# DEV by Shuhuai Li

> [!important]
> Under development. Just for reference

The script under `./hotspot` gets the time-consuming phases thoughout the workflow. Now, I focus on the `hisat-3n` tool. You can check my modification on codes [here](https://github.com/ASC25-SHU02/hisat2/tree/getline_op). 
I name my work after branch name `XXX_op`, where XXX explains what aspect I try to optimize in this branch, while op stands for optimization (genshin impact).

## TODO
1. Overall workflow optimization:
    - four hotspot commands run in parallel
    - substitue `unix pipe` for temporary files
2. Lock-free queue
在C++中，实现单生产者多消费者的无锁队列有多个成熟的库可供选择。以下是几个常用的库：

### 1. **Boost.Lockfree**
   - **简介**: Boost库提供了一个无锁队列的实现，支持单生产者单消费者（SPSC）和单生产者多消费者（SPMC）模式。
   - **特点**: 高性能，易于集成，适合需要无锁队列的场景。
   - **示例**:
     ```cpp
     #include <boost/lockfree/queue.hpp>
     #include <thread>
     #include <iostream>

     boost::lockfree::queue<int> queue(128);

     void producer() {
         for (int i = 0; i < 100; ++i) {
             while (!queue.push(i)) { }
         }
     }

     void consumer(int id) {
         int value;
         while (queue.pop(value)) {
             std::cout << "Consumer " << id << " popped " << value << std::endl;
         }
     }

     int main() {
         std::thread prod(producer);
         std::thread cons1(consumer, 1);
         std::thread cons2(consumer, 2);

         prod.join();
         cons1.join();
         cons2.join();

         return 0;
     }
     ```

### 2. **Moodycamel::ConcurrentQueue**
   - **简介**: 这是一个高性能的无锁队列，支持多生产者多消费者（MPMC）模式，同时也适用于单生产者多消费者（SPMC）场景。
   - **特点**: 高性能，易于使用，支持批量操作。
   - **示例**:
     ```cpp
     #include "concurrentqueue.h"
     #include <thread>
     #include <iostream>

     moodycamel::ConcurrentQueue<int> queue;

     void producer() {
         for (int i = 0; i < 100; ++i) {
             queue.enqueue(i);
         }
     }

     void consumer(int id) {
         int value;
         while (queue.try_dequeue(value)) {
             std::cout << "Consumer " << id << " popped " << value << std::endl;
         }
     }

     int main() {
         std::thread prod(producer);
         std::thread cons1(consumer, 1);
         std::thread cons2(consumer, 2);

         prod.join();
         cons1.join();
         cons2.join();

         return 0;
     }
     ```

### 3. **Folly's MPMCQueue**
   - **简介**: Folly是Facebook开源的C++库，其中的`MPMCQueue`是一个高性能的无锁队列，支持多生产者多消费者模式。
   - **特点**: 高性能，适合高并发场景。
   - **示例**:
     ```cpp
     #include <folly/MPMCQueue.h>
     #include <thread>
     #include <iostream>

     folly::MPMCQueue<int> queue(128);

     void producer() {
         for (int i = 0; i < 100; ++i) {
             queue.write(i);
         }
     }

     void consumer(int id) {
         int value;
         while (queue.read(value)) {
             std::cout << "Consumer " << id << " popped " << value << std::endl;
         }
     }

     int main() {
         std::thread prod(producer);
         std::thread cons1(consumer, 1);
         std::thread cons2(consumer, 2);

         prod.join();
         cons1.join();
         cons2.join();

         return 0;
     }
     ```

### 4. **TBB (Intel Threading Building Blocks)**
   - **简介**: Intel TBB库提供了一个并发队列`concurrent_queue`，支持多生产者多消费者模式。
   - **特点**: 高性能，适合多线程编程。
   - **示例**:
     ```cpp
     #include <tbb/concurrent_queue.h>
     #include <thread>
     #include <iostream>

     tbb::concurrent_queue<int> queue;

     void producer() {
         for (int i = 0; i < 100; ++i) {
             queue.push(i);
         }
     }

     void consumer(int id) {
         int value;
         while (queue.try_pop(value)) {
             std::cout << "Consumer " << id << " popped " << value << std::endl;
         }
     }

     int main() {
         std::thread prod(producer);
         std::thread cons1(consumer, 1);
         std::thread cons2(consumer, 2);

         prod.join();
         cons1.join();
         cons2.join();

         return 0;
     }
     ```

### 总结
- **Boost.Lockfree**: 适合需要轻量级无锁队列的场景。
- **Moodycamel::ConcurrentQueue**: 高性能，支持批量操作，适合高并发场景。
- **Folly's MPMCQueue**: 高性能，适合高并发场景，尤其是Facebook的开源项目。
- **TBB**: 适合需要与Intel TBB其他组件集成的场景。

---
# Origin stuff (Haven't changed)

# m<sup>5</sup>C-UBSseq

## Changelog

- 4/23/2024: rewrite code using polars

## workflow

[![](./docs/flow.svg)](https://github.com/y9c/m5C-UBSseq)

## Citation

- cite this software

  ```BibTex
  @software{chang_y_2024_11046885,
      author    = {Chang Y},
      title     = {y9c/m5C-UBSseq: V0.1},
      publisher = {Zenodo},
      version   = {v0.1},
      doi       = {10.5281/zenodo.11046885},
      url       = {https://doi.org/10.5281/zenodo.11046885}
  }
  ```

- cite the method

  ```BibTex
  @article{dai_ultrafast_2024,
      title = {Ultrafast bisulfite sequencing detection of 5-methylcytosine in {DNA} and {RNA}},
      url = {https://www.nature.com/articles/s41587-023-02034-w},
      doi = {10.1038/s41587-023-02034-w},
      author = {Dai, Qing and Ye, Chang and Irkliyenko, Iryna and Wang, Yiding and Sun, Hui-Lung and Gao, Yun and Liu, Yushuai and Beadell, Alana and Perea, José and Goel, Ajay and He, Chuan},
      date = {2024-01-02},
  }
  ```

&nbsp;

<p align="center">
<img
  src="https://raw.githubusercontent.com/y9c/y9c/master/resource/footer_line.svg?sanitize=true"
/>
</p>
<p align="center">
Copyright &copy; 2021-present
<a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
