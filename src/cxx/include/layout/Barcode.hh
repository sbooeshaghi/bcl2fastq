/**
 * BCL to FASTQ file converter
 * Copyright (c) 2007-2017 Illumina, Inc.
 *
 * This software is covered by the accompanying EULA
 * and certain third party copyright/licenses, and any user of this
 * source file is bound by the terms therein.
 *
 * \file Barcode.hh
 *
 * \brief Declaration of barcode.
 *
 * \author Marek Balint
 */


#ifndef BCL2FASTQ_LAYOUT_BARCODE_HH
#define BCL2FASTQ_LAYOUT_BARCODE_HH


#include <vector>
#include <string>

#include <stdint.h>

#include "common/Types.hh"

namespace bcl2fastq {


namespace layout {


/// \brief Barcode type definition.
class Barcode
{
public:

    /// \brief Barcode component type definition.
    class Component
    {
    public:

        /// \brief Size type definition.
        /// \note Size of the integer type limits maximum number of component's bases (1 base per 3 bits).
        typedef uint8_t SizeType;

    public:

        /// \brief Invalid base character.
        static const char INVALID_BASE;

    public:

        /// \brief Constructor.
        /// \param value Component value.
        /// \pre <code>value.empty() == false</code>
        /// \pre <code>value.size() <= 21</code>
        explicit Component(const std::string &value);

        Component() : value_() { }

    public:

        /// \brief Less-than comparison.
        /// \param rhs Other instance to compare with.
        /// \retval true Other instance is less than this instance.
        /// \retval false Other instance is not less than this instance.
        bool less(const Component &rhs) const { return value_ < rhs.value_; }

        bool operator!=(const Component &rhs) const { return value_ != rhs.value_; }
    public:

        /// \brief Get barcode component length.
        /// \return Number of bases in barcode component.
        Barcode::Component::SizeType getLength() const { return value_.size(); }

        /// \brief Get the barcode component length.
        /// \return Number of bases in barcode component.
        size_t size() const { return value_.size(); }

        /// \brief Get the barcode string.
        /// \return the string.
        const std::string& getString() const { return value_; }

        /// \brief Get the barcode string.
        /// \return the string.
        std::string& getString() { return value_; }

        /// \brief Reset the memory, without deallocating.
        void reset() { value_.clear(); }

    public:

        /// \brief Generate components with mismatches.
        /// \param position Mismatching position.
        /// \param componentInserter Insert iterator for mismatching barcode components.
        /// \pre <code>position < this->getLength()</code>
        template<typename Inserter>
        void generateMismatches(
            Component::SizeType position,
            Inserter componentInserter
        ) const;

    private:

        /// \brief Value.
        std::string value_;
    };

    /// \brief Barcodes container type definition.
    typedef std::vector<Component> ComponentsContainer;

public:

    /// \brief Barcode components separator.
    static const char COMPONENT_SEPARATOR;

    /// \brief Label to use when no Barcode is given.
    static const char DEFAULT_BARCODE[];

public:

    /// \brief Constructor.
    /// \param componentsBegin Beginning of barcode components.
    /// \param componentsEnd End of barcode components.
    /// \pre Iterator @c componentsEnd is reachable from @c componentsBegin.
    /// \pre <code>componentsBegin != componentsEnd</code>
    template<typename Iterator>
    Barcode(Iterator componentsBegin, Iterator componentsEnd);

    Barcode() : components_() { } 

    void reset(size_t numComponents);

    std::string toString() const;

public:

    /// \brief Less-than comparison.
    /// \param rhs Other instance to compare with.
    /// \retval true Other instance is less than this instance.
    /// \retval false Other instance is not less than this instance.
    bool operator<(const Barcode &rhs) const;

    bool operator==(const Barcode &rhs) const;
    bool operator!=(const Barcode &rhs) const;

public:

    /// \brief Get beginning of barcode components.
    /// \return Iterator to beginning of barcode segmetns.
    Barcode::ComponentsContainer::const_iterator componentsBegin() const { return components_.begin(); }

    /// \brief Get end of barcode components.
    /// \return Iterator to end of barcode components.
    Barcode::ComponentsContainer::const_iterator componentsEnd() const { return components_.end(); }

    /// \brief Get the barcode components
    /// \return Barcode components
    const ComponentsContainer& getComponents() const { return components_; }

public:

    /// \brief Generate barcodes with mismatches.
    /// \param component Mismatching component.
    /// \param position Mismatching position of component.
    /// \param barcodeInserter Insert iterator for mismatching barcodes.
    /// \pre Iterator @c component must be reachable from <code>this->componentsBegin()</code>
    /// \pre Iterator @c component must be dereferencable.
    /// \pre <code>position < component->getLength()</code>
    template<typename Inserter>
    void generateMismatches(
        ComponentsContainer::const_iterator component,
        Component::SizeType position,
        Inserter barcodeInserter
    ) const;

    /// \brief Ouput the barcode to the buffer
    /// \param buffer to output to
    void output(std::vector<char>& buffer) const;

    /// \brief Mask the barcode for the given index number
    /// \param index number.
    void mask(common::ReadNumber indexNumber);

//private:

    friend std::ostream& operator<<(std::ostream& os, const Barcode& barcode);

    /// \brief Barcode components.
    ComponentsContainer components_;
};

/// \brief Barcodes container type definition.
typedef std::vector<Barcode> BarcodesContainer;

} // namespace layout


} // namespace bcl2fastq


#include "layout/Barcode.hpp"


#endif // BCL2FASTQ_LAYOUT_BARCODE_HH


