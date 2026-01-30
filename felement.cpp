#include "felement.h"

FElement::FElement(QObject *parent)
    : QObject{parent}
{}

const unsigned short &FElement::getElementNumber() const
{
    return elementNumber;
}

void FElement::setElementNumber(unsigned short newElementNumber)
{
    elementNumber = newElementNumber;
}

const QString &FElement::getElementName() const
{
    return elementName;
}

void FElement::setElementName(const QString &newElementName)
{
    elementName = newElementName;
}

const QColor &FElement::getElementColor() const
{
    return elementColor;
}

void FElement::setElementColor(const QColor &newElementColor)
{
    elementColor = newElementColor;
}
