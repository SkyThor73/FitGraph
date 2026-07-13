#ifndef FDATAPCS_H
#define FDATAPCS_H

#include "foutput.h"
#include "felement.h"

class FDataPCS
{
public:
    FDataPCS();
    static const QString pcs_path;                              //Путь файла с сечениями
    static void initializeDatabase();                           //Глобальный метод инициализации базы данных
    static unsigned short getAtomicNumber(const QString &symbol);//Глобальный метод извлечения зарядового числа
    enum Shells { S = 'S', P = 'P', D = 'D', F = 'F', G = 'G' }; //Перечисление оболочек
    void MiningPCS(const FOutput &output, QList<FElement *> &element);

    const QList<QList<QList<QList<double> > > > &getPhotoionizationCrossSection() const;

    const QList<QList<QList<QList<double> > > > &getAsymmetry() const;

private:
    QList<QList<QList<QList<double>>>> photoionizationCrossSection, asymmetry;
};

#endif // FDATAPCS_H
