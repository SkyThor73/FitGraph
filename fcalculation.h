#ifndef FCALCULATION_H
#define FCALCULATION_H

#include <QObject>
#include "foutput.h"
#include "fdatapcs.h"

class FCalculation : public QObject
{
    Q_OBJECT
public:
    explicit FCalculation(QObject *parent = nullptr);
    enum Sources { HE_I, HE_II, MG, AL };   //Перечисление источников
    void Calculation(const FOutput &output, const QList<FElement *> &element);    //Выполняет расчёт и заполнение всех spd векторов
    void CalculationIntegralArea(const FOutput &output, const FDataPCS &pcs, const unsigned short &set_orbitals, const QList<double> &setEnergy, const QList<QList<QList<QList<double> > > > &spdSet);  //Вызывается при изменении source и resolution_ratio
    void CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy);  //Вызывается при изменении number_of_steps, initial_value и end_value
    void CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy, QList<double> resolution_vector);  //Вызывается при изменении number_of_steps, initial_value и end_value
    void CalculationEachWidth(const FOutput &output, const FDataPCS &pcs, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &spdSet, QList<double> resolution_vector);

    const Sources &getSource() const;
    void setSource(Sources newSource);

    const unsigned short &getSteps() const;
    void setSteps(unsigned short newSteps);

    const double &getInitial_value() const;
    void setInitial_value(double newInitial_value);

    const double &getEnd_value() const;
    void setEnd_value(double newEnd_value);

    const double &getResolution_ratio() const;
    void setResolution_ratio(double newResolution_ratio);

    const double &getMax_value_y() const;

    const double &getUpper_bound() const;
    void setUpper_bound(double newUpper_bound);

    const QList<double> &getIntegral_area() const;

    const QList<double> &getAbscissaX() const;

    const QList<double> &getOrdinateY() const;

    const QList<double> &getResolution_vector() const;
    void setResolution_vector(const QList<double> &newResolution_vector);

    const QList<double> &getRIntensity() const;

    const QList<QList<double> > &getGaussian() const;

    const QList<QList<QList<QList<double> > > > &getSpdSimple() const;

    const QList<QList<QList<QList<double> > > > &getSpdSimpleBeta() const;

    const QList<QList<QList<QList<double> > > > &getSpdMulliken() const;

    const QList<QList<QList<QList<double> > > > &getSpdMullikenBeta() const;

    const QList<QList<QList<QList<double> > > > &getSpdLowdin() const;

    const QList<QList<QList<QList<double> > > > &getSpdLowdinBeta() const;

    const QList<QList<QList<QList<double> > > > &getSpdRIntensity() const;

signals:

private:
    Sources source = HE_I;  //Переменная источника
    unsigned short steps = 0;
    double initial_value = 0.0, end_value = 0.0, resolution_ratio = 0.01, max_value_y = 0.0, upper_bound = 0.0;
    QList<double> integral_area, abscissaX, ordinateY, resolution_vector, RIntensity, betaMO;
    QList<QList<double> > gaussian;
    QList<QList<QList<QList<double> > > > spdSimple, spdSimpleBeta, spdMulliken, spdMullikenBeta, spdLowdin, spdLowdinBeta, spdRIntensity;

    void RefreshVector(const FOutput &output, const unsigned short &set_orbitals, QList<QList<QList<QList<double> > > > &spd);  //Обновление вектора при выборе нового файла
    void FillingVector(const FOutput &output, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &uVectors, const QList<double> &sumVector, QList<QList<QList<QList<double> > > > &spdVectorAB); //Заполнение вектора
    QList<double> sumSPD(const FOutput &output, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &uVectors);
    QList<QList<QList<QList<double> > > > ConvolusionSPD(const FOutput &output, const QList<FElement *> &element, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &setVectors);
};

#endif // FCALCULATION_H
