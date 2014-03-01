/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class HasClassRegistry -- implementation of HasClassRegistry template class.
*   This file should be included from just one cpp file and instantiated using
*   template class HasClassRegistry<ConcreteBase>;
*   This ensure there is just one instance of the RegistryStorage
*
*****************************************************************************/

#ifndef HASCLASSREGISTRY_IPP_INCLUDED
#define HASCLASSREGISTRY_IPP_INCLUDED

#include <memory>

namespace diffpy {

// Public Methods ------------------------------------------------------------

template <class TBase>
bool HasClassRegistry<TBase>::registerThisType() const
{
    using namespace std;
    RegistryStorage& reg = getRegistry();
    if (reg.count(this->type()))
    {
        // do nothing when registering the same class twice
        const TBase& regprot = *(reg[this->type()]);
        if (typeid(*this) == typeid(regprot))    return true;
        // raise exception if trying to register a different class
        ostringstream emsg;
        emsg << "Prototype type '" << this->type() <<
            "' is already registered.";
        throw logic_error(emsg.str());
    }
    SharedPtr p = this->create();
    this->setupRegisteredObject(p);
    reg[this->type()] = p;
    return true;
}

// Public Static Methods -----------------------------------------------------

template <class TBase>
bool HasClassRegistry<TBase>::aliasType(const std::string& tp,
        const std::string& al)
{
    using namespace std;
    RegistryStorage& reg = getRegistry();
    if (!reg.count(tp))
    {
        ostringstream emsg;
        emsg << "Cannot create alias for unknown prototype '" <<
            tp << "'.";
        throw logic_error(emsg.str());
    }
    if (reg.count(al) && reg[al] != reg[tp])
    {
        ostringstream emsg;
        emsg << "Prototype type '" << al <<
            "' is already registered.";
        throw logic_error(emsg.str());
    }
    reg[al] = reg[tp];
    return true;
}


template <class TBase>
typename HasClassRegistry<TBase>::SharedPtr
HasClassRegistry<TBase>::createByType(const std::string& tp)
{
    using namespace std;
    typename RegistryStorage::iterator irg;
    RegistryStorage& reg = getRegistry();
    irg = reg.find(tp);
    if (irg == reg.end())
    {
        ostringstream emsg;
        emsg << "Unknown type '" << tp << "'.";
        throw invalid_argument(emsg.str());
    }
    SharedPtr rv = irg->second->create();
    return rv;
}


template <class TBase>
std::set<std::string>
HasClassRegistry<TBase>::getRegisteredTypes()
{
    using namespace std;
    set<string> rv;
    RegistryStorage& reg = getRegistry();
    typename RegistryStorage::iterator irg;
    for (irg = reg.begin(); irg != reg.end(); ++irg)
    {
        rv.insert(irg->second->type());
    }
    return rv;
}

// Private Static Methods ----------------------------------------------------

template <class TBase>
typename HasClassRegistry<TBase>::RegistryStorage&
HasClassRegistry<TBase>::getRegistry()
{
    static std::auto_ptr<RegistryStorage> the_registry;
    if (!the_registry.get())
    {
        the_registry.reset(new RegistryStorage());
    }
    return *the_registry;
}

}   // namespace diffpy

// Make these definitions easily accessible from including cpp file.

using diffpy::HasClassRegistry;

// vim:ft=cpp:

#endif  // HASCLASSREGISTRY_IPP_INCLUDED
