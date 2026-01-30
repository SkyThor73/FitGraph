#ifndef FEXPERIMENTAL_H
#define FEXPERIMENTAL_H

#include <QObject>
#include "fcalculation.h"

class FExperimental : public QObject
{
    Q_OBJECT
public:
    explicit FExperimental(QObject *parent = nullptr);
    void ReadXY(const FCalculation &calculate);


    const double &getMax_value_exp_y() const;
    void setMax_value_exp_y(double newMax_value_exp_y);

    const QString &getPath_exp_file() const;
    void setPath_exp_file(const QString &newPath_exp_file);

    const QList<double> &getExpOrdinateY() const;
    void setExpOrdinateY(const QList<double> &newExpOrdinateY);

    const QList<double> &getExpAbscissaX() const;

signals:

private:
    double max_value_exp_y = 0.0;
    QString path_exp_file = "";
    QList<double> expAbscissaX, expOrdinateY;
};

#endif // FEXPERIMENTAL_H
