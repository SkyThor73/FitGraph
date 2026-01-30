#ifndef FELEMENT_H
#define FELEMENT_H

#include "qcolor.h"
#include <QObject>

class FElement : public QObject
{
    Q_OBJECT
public:
    explicit FElement(QObject *parent = nullptr);

    const unsigned short &getElementNumber() const;
    void setElementNumber(unsigned short newElementNumber);

    const QString &getElementName() const;
    void setElementName(const QString &newElementName);

    const QColor &getElementColor() const;
    void setElementColor(const QColor &newElementColor);

private:
    unsigned short elementNumber;
    QString elementName;
    QColor elementColor;

signals:
};

#endif // FELEMENT_H
