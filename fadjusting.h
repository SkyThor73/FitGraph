#ifndef FADJUSTING_H
#define FADJUSTING_H

#include <QObject>
#include "fcalculation.h"
#include "fexperimental.h"

class FAdjusting : public QObject
{
    Q_OBJECT
public:
    explicit FAdjusting(QObject *parent = nullptr);
    void FindBetterShift(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp);
    void FindBetterVise(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp);
    void FindBetterWidth(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp);
    void ShiftEnergy(QList<double> &new_energy);
    void CalculationCompliance(const FCalculation &calculate, const FExperimental &exp);    //Расчёт соответствия

    const double &getSpecific_min_square() const;

    const double &getVisible_shift_x() const;
    void setVisible_shift_x(double newVisible_shift_x);

signals:

private:
    double shift_x = 0.0, visible_shift_x = 0.0, specific_min_square = 0.0, fwhm_min = 0.01, fwhm_max = 2.5;
};

#endif // FADJUSTING_H
