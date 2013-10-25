/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class Attributes - interface for calling setter and getter methods using
*   their string names.
*
*****************************************************************************/

#ifndef ATTRIBUTES_HPP_INCLUDED
#define ATTRIBUTES_HPP_INCLUDED

#include <string>
#include <set>
#include <map>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

namespace diffpy {
namespace attributes {

class Attributes;

/// @class DoubleAttributeError
/// @brief custom exception for Attributes-related errors.  This is
/// thrown for invalid names or for attempts to set a read-only attribute.

class DoubleAttributeError : public std::runtime_error
{
    public:
        DoubleAttributeError(const std::string msg="") :
            std::runtime_error(msg)
        { }
};

/// @class BaseDoubleAttribute
/// @brief abstract base class for accessing a particular double attribute

class BaseDoubleAttribute
{
    public:

        virtual ~BaseDoubleAttribute() { }
        virtual double getValue(const Attributes* obj) const = 0;
        virtual void setValue(Attributes* obj, double value) = 0;
        virtual bool isreadonly() const = 0;
};


class BaseAttributesVisitor
{
    public:

        virtual ~BaseAttributesVisitor() { }

        virtual void visit(const Attributes& a) = 0;
        virtual void visit(Attributes& a)
        {
            visit(static_cast<const Attributes&>(a));
        }
};

/// @class Attributes
/// @brief implementation of attribute access.  The client classes
/// should derive from Attributes and register their setter and
/// getter methods in their constructors.

class Attributes
{
    public:

        // class is virtual
        virtual ~Attributes()  { }

        // assignments in derived classes should not change mdoubleattrs
        // mdoubleattrs should be only changed via registerDoubleAttribute
        Attributes& operator=(const Attributes& other)  { return *this; }

        // methods
        double getDoubleAttr(const std::string& name) const;
        void setDoubleAttr(const std::string& name, double value);
        bool hasDoubleAttr(const std::string& name) const;
        std::set<std::string> namesOfDoubleAttributes() const;
        std::set<std::string> namesOfWritableDoubleAttributes() const;
        // visitors
        virtual void accept(BaseAttributesVisitor& v)  { v.visit(*this); }
        virtual void accept(BaseAttributesVisitor& v) const  { v.visit(*this); }

    protected:

        friend void registerBaseDoubleAttribute(Attributes*,
                const std::string&, attributes::BaseDoubleAttribute* pa);
        template <class T, class Getter>
            void registerDoubleAttribute(const std::string& name, T* obj, Getter);
        template <class T, class Getter, class Setter>
            void registerDoubleAttribute(const std::string& name, T* obj, Getter, Setter);

    private:

        // types
        typedef std::map<std::string,
                boost::shared_ptr<attributes::BaseDoubleAttribute> >
                    DoubleAttributeStorage;
        // data
        DoubleAttributeStorage mdoubleattrs;

        // methods
        void checkAttributeName(const std::string& name) const;

        // visitor classes

        class CountDoubleAttrVisitor : public BaseAttributesVisitor
        {
            public:

                CountDoubleAttrVisitor(const std::string& name);
                virtual void visit(const Attributes& a);
                int count() const;

            private:

                // data
                const std::string& mname;
                int mcount;
        };


        class GetDoubleAttrVisitor : public BaseAttributesVisitor
        {
            public:

                GetDoubleAttrVisitor(const std::string& name);
                virtual void visit(const Attributes& a);
                double getValue() const;

            private:

                // data
                const std::string& mname;
                double mvalue;
        };


        class SetDoubleAttrVisitor : public BaseAttributesVisitor
        {
            public:

                SetDoubleAttrVisitor(const std::string& name, double value);
                virtual void visit(const Attributes& a);
                virtual void visit(Attributes& a);

            private:

                // data
                const std::string& mname;
                double mvalue;
        };


        class NamesOfDoubleAttributesVisitor : public BaseAttributesVisitor
        {
            public:

                NamesOfDoubleAttributesVisitor(bool excludereadonly);
                virtual void visit(const Attributes& a);
                const std::set<std::string>& names() const;

            private:

                // data
                bool mexcludereadonly;
                std::set<std::string> mnames;
        };

};  // class Attributes

// non-member helpers

void registerBaseDoubleAttribute(Attributes* obj,
        const std::string& name, BaseDoubleAttribute* pa);

void throwDoubleAttributeReadOnly();

typedef std::map<std::string, double> AttributesDataMap;

AttributesDataMap saveAttributesData(const Attributes& obj);

void loadAttributesData(Attributes& obj, const AttributesDataMap& data);

}   // namespace attributes
}   // namespace diffpy

// Implementation ------------------------------------------------------------

#include <diffpy/Attributes.ipp>

// make selected classes visible in diffpy namespace
namespace diffpy {
    using attributes::Attributes;
    using attributes::BaseAttributesVisitor;
}

#endif  // ATTRIBUTES_HPP_INCLUDED
