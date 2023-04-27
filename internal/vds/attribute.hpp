#ifndef VDS_SLICE_ATTRIBUTE_HPP
#define VDS_SLICE_ATTRIBUTE_HPP

#include <functional>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include "regularsurface.hpp"

template< typename Derived >
class Attribute {
public:
    Attribute(void* dst, std::size_t size) : dst(dst), size(size) {};

    template< typename InputIt >
    void compute(InputIt begin, InputIt end) noexcept (false) {
        static_cast< Derived& >(*this).compute(begin, end);
    }

    void fill(float x) noexcept (true) { this->write(x); };

private:
    void*       dst;
    std::size_t size;
    std::size_t offset = 0;

protected:
    void write(float x) {
        char* dst = (char*)this->dst + this->offset * sizeof(float);
        memcpy(dst, &x, sizeof(float));
        ++this->offset;
    }
};

struct Min : public Attribute< Min > {
    Min(void* dst, std::size_t size) : Attribute< Min >(dst, size) {}

    template< typename InputIt >
    void compute(InputIt begin, InputIt end) noexcept (false) {
        float value =  *std::min_element(begin, end);
        this->write(value);
    }
};

struct Max : public Attribute< Max > {
    Max(void* dst, std::size_t size) : Attribute< Max >(dst, size) {}

    template< typename InputIt >
    void compute(InputIt begin, InputIt end) noexcept (false) {
        float value = *std::max_element(begin, end);
        this->write(value);
    }
};

struct Mean : public Attribute< Mean > {
    Mean(void* dst, std::size_t size, std::size_t vsize)
        : Attribute< Mean >(dst, size)
        , vsize(vsize)
    {}

    template< typename InputIt >
    void compute(InputIt begin, InputIt end) noexcept (false) {
        float sum = std::accumulate(begin, end, 0.0f);
        this->write(sum / this->vsize);
    }
private:
    std::size_t vsize;
};

struct Rms : public Attribute< Rms > {
    Rms(void* dst, std::size_t size, std::size_t vsize)
        : Attribute< Rms >(dst, size)
        , vsize(vsize)
    {}

    template< typename InputIt >
    void compute(InputIt begin, InputIt end) noexcept (false) {
        float sum = std::accumulate( std::next(begin), end, std::pow(*begin, 2),
            [](float a, float b) {
                return a + std::pow(b, 2);
            }
        );
        this->write(std::sqrt(sum / this->vsize));
    }
private:
    std::size_t vsize;
};

using Attributes = std::variant< Min, Max, Mean, Rms >;

struct AttributeFillVisitor {
    AttributeFillVisitor( float fillvalue ) : fillvalue(fillvalue) {}

    template< typename Attr >
    void operator()(Attr& attribute) const {
        attribute.fill( this->fillvalue );
    }

private:
    float fillvalue;
};

template< typename It >
struct AttributeComputeVisitor {
    AttributeComputeVisitor(It begin, It end) : begin(begin), end(end) {}

    template< typename Attr >
    void operator()(Attr& attribute) const {
        attribute.compute(this->begin, this->end);
    }

private:
    It begin;
    It end;
};


/**
 * A Window is defined as some area around a reference point.
 */
struct Window {
    Window(float up, float down, float samplerate)
        : m_up(up), m_down(down), m_samplerate(samplerate)
    {
        // Verify that consistency
    }

    float distance_up() const { return this->m_up; }
    float distance_down() const { return this->m_down; }

    float samples_up()   const { return this->m_up   / this->m_samplerate; }
    float samples_down() const { return this->m_down / this->m_samplerate; }

    float samplerate() const { return this->m_samplerate; }
    std::size_t size() const { return this->samples_up() + 1 + this->samples_down(); }

private:
    float m_samplerate;
    float m_up;
    float m_down;
};

template< typename Derived >
struct WindowResampler {

    WindowResampler(Window input, Window output)
        : input( std::move(input) )
        , output( std::move(output) )
    {}

    std::size_t srcsize() const { return this->input.nsamples(); }
    std::size_t dstsize() const { return this->output.nsamples(); }

    template< typename OutputIt >
    void resample(void const * src, std::size_t srcsize, OutputIt begin, OutputIt end) const {
        static_cast< Derived >(*this)->resample(src, srcsize, begin, end);
    }
protected:
    Window input;
    Window output;
};

struct CubicWindowResampler : public WindowResampler< CubicWindowResampler > {
    CubicWindowResampler(Window input, Window output)
        : WindowResampler< CubicWindowResampler >(input, output)
    {}

    template< typename OutputIt >
    void resample(
        void const* src,
        std::size_t srcsize,
        OutputIt begin,
        OutputIt end,
        float offset
    ) const {
        boost::math::interpolators::cardinal_cubic_b_spline< float > spline(
            static_cast< float const* >( src ),
            this->input.size(),
            0,
            this->input.samplerate()
        );
        float i = offset;
        std::transform(begin, end, begin,
            [&](float const& n){
                auto v = spline(i);
                i += this->output.samplerate();
                return v;
            }
        );
    }
};

/** Windowed horizon
 *
 * Definition and layout of a windowed horizon
 * -------------------------------------------
 *
 * A windowed horizon is a horizon where we have multiple vertical
 * sample values for each horizontal position. Typically defined by a window
 * above and a window below the horizon definition. As illustrated by this
 * *flat* horizon:
 *                       *------*
 *                      /      /|
 *                     /      / |
 *                  - *------* /|
 *     window above | |      |/ *
 * horizon ---->    - |------| /
 *     window below | |      |/
 *                  - *------*
 *
 *
 * Note that the above and below window is just and example of how to think
 * about this data structure. This class, however, doesn't care about the
 * window definition in that much detail. It only needs to know the size of the
 * full vertical window. I.e. how many samples there are at every horizontal
 * position.
 *
 * As the illustration show, the windowed horizon is a 3D array. It's assumed
 * that the vertical axis is the fastest moving. I.e. that vertical samples at
 * the same horizontal position are contiguous in memory. There are no
 * assumptions on the order of the 2 horizontal axes.
 *
 * Class features
 * --------------
 *
 * In short this class is an abstraction around a float pointer to a 3D array
 * with the properties described above. The reason for the raw pointer is to be
 * compatible with C, and hence CGO and GO.
 *
 * `calc_attribute` is the main interface and provides a way for the caller to
 * specify a calculation (in form of a lambda) that is applied to every
 * vertical window. This makes it very easy for the consumer to implement their
 * own attributes without having to implement correct looping of the windows
 * each time.
 */
class Horizon{
private:
    struct StridedIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = float;
        using pointer           = const float*;
        using reference         = const float&;

        StridedIterator(pointer cur, std::size_t step) : cur(cur), step(step) {}

        reference operator*()  { return *this->cur; }
        pointer   operator->() { return  this->cur; }

        StridedIterator& operator++() { this->cur += this->step; return *this; };

        friend
        bool operator==(StridedIterator const& lhs, StridedIterator const& rhs) {
            return lhs.cur == rhs.cur and lhs.step == rhs.step;
        }

        friend
        bool operator!=(StridedIterator const& lhs, StridedIterator const& rhs) {
            return !(lhs == rhs);
        }
    private:
        pointer cur;
        std::size_t step;
    };

    struct VerticalIterator : public StridedIterator {
        VerticalIterator(StridedIterator::pointer cur) : StridedIterator(cur, 1) {}
    };

public:

    Horizon(
        const float*   data,
        RegularSurface surface,
        Window         vertical,
        float          fillvalue
    ) : m_ptr(data)
      , m_surface( std::move(surface) )
      , m_vertical( std::move(vertical) )
      , m_fillvalue(fillvalue)
    {}

    using HorizontalIt = StridedIterator;
    using VerticalIt   = VerticalIterator;

    HorizontalIt begin() const noexcept (true);
    HorizontalIt end()   const noexcept (true);

    RegularSurface surface() const noexcept (true) { return this->m_surface; };
    Window vertical_window() const noexcept (true) { return this->m_vertical; };

    std::size_t size() const noexcept (true) {
        return this->m_surface.size() * this->m_vertical.size();
    }

    float fillvalue() const noexcept (true) { return this->m_fillvalue; };

    /** Workhorse for attribute calculation
     *
     * The provided lambda func is called with an iterator pair (begin, end)
     * for each vertical window in the horizon. The results of each lambda call
     * is memcpy'd into the output-buffer 'dst'. 'dst' needs to be at least
     * this->mapsize() large. Example of an lambda that will calculate the
     * minimum value of each vertical window:
     *
     *     auto min_value = [](It beg, It end) {
     *         return *std::min_element(beg, end);
     *     }
     *
     * Calling calc_attribute() like so, will write the entire attribute to dst:
     *
     *     horizon.calc_attribute(dst, size, min_value)
     *
     * Like the horzion itself, attribute values are written as floating
     * points.
     *
     * The first value of each window is compared to fillvalue. If true, no
     * computation takes place for that window and the output buffer will be
     * populated with the fillvalue instead.
     */
    template< typename Func >
    void calc_attribute(void* dst, std::size_t size, Func func) const;

    template< typename Resampler >
    void calc_attributes(
        HorizontalIt begin,
        HorizontalIt end,
        Window target_window,
        std::vector< Attributes > & attrs
    ) const noexcept (false) {
        std::vector< float > buffer(target_window.size());

        Resampler resampler(this->vertical_window(), target_window);

        AttributeFillVisitor    fillvisitor(this->fillvalue());
        AttributeComputeVisitor computevisitor(buffer.begin(), buffer.end());
        //TODO assert > 0. Or have an horizon.window.encapulate(target_window);

        auto i = std::distance(this->begin(), begin);
        auto compute = [&](const float& front) {
            if (front == this->fillvalue()) {
                std::for_each(attrs.begin(), attrs.end(), [&](Attributes& attr) {
                    std::visit(fillvisitor, attr);
                });
            } else {

                float ref = this->m_surface.at(i);
                float window_offset = std::fmod(ref, this->vertical_window().samplerate());

                resampler.resample(
                    &front,
                    this->m_vertical.size(),
                    buffer.begin(),
                    buffer.end(),
                    window_offset
                );

                std::for_each(attrs.begin(), attrs.end(), [&](Attributes& attr) {
                    std::visit(computevisitor, attr);
                });
            }

            ++i;
        };

        std::for_each(begin, end, compute);
    }

private:
    const          float* m_ptr;
    RegularSurface m_surface;
    Window         m_vertical;
    float          m_fillvalue;
};

#endif /* VDS_SLICE_ATTRIBUTE_HPP */
