#ifndef FOUTPUT_H
#define FOUTPUT_H

#include <QObject>

class FOutput : public QObject
{
    Q_OBJECT
public:
    explicit FOutput(QObject *parent = nullptr);
    enum Programms { NONE, GAMESS, FIREFLY, ORCA, GAUSSIAN };               //Перечисление программ квантовой химии
    void ReadData(const QString &path_file);                                //Считывает все данные из файла и присваивает их в объект класса
    void ReadingOrca(QTextStream &stream, QString &str);
    void ReadingGamessFly(QTextStream &stream, QString &str);
    void ReadingGaussian(QTextStream &stream, QString &str);

    const Programms &getProgramm() const;

    const bool &getUhf() const;

    const bool &getMulliken() const;

    const bool &getLowdin() const;

    const unsigned short &getKSpdfghi() const;

    const unsigned short &getOrbital_count() const;

    const unsigned short &getOrbital_count_beta() const;

    const unsigned short &getAtom_count() const;

    const unsigned short &getUnique_atoms() const;

    const QString &getBasis() const;

    const QList<double> &getEnergy() const;

    const QList<double> &getEnergyBeta() const;

    const QList<QString> &getElements() const;

    const QList<QString> &getUElements() const;

    const QList<QString> &getSymmetry() const;

    const QList<QString> &getSymmetryBeta() const;

    const QList<QList<QList<QList<double> > > > &getEigenVectors() const;

    const QList<QList<QList<QList<double> > > > &getEigenVectorsBeta() const;

    const QList<QList<QList<QList<double> > > > &getMullikenVectors() const;

    const QList<QList<QList<QList<double> > > > &getMullikenVectorsBeta() const;

    const QList<QList<QList<QList<double> > > > &getLowdinVector() const;

    const QList<QList<QList<QList<double> > > > &getLowdinVectorBeta() const;

signals:

private:
    const double kHartree = 27.2113862459885353;        //Постоянная хартри
    const unsigned short kSpdfghi = 6;                  //Количество атомных подуровней

    Programms programm = NONE;

    bool runtype = false, mulliken = false, uhf = false, lowdin = false, ecp = false, scftyp = false, cctyp = false;
    unsigned short basis_functions = 0, orbital_count = 0, orbital_count_beta = 0, atom_count = 0, unique_atoms = 0;
    QString basis = "";
    QList<bool> ecpAtom;
    QList<unsigned short> ecpElectrons;
    QList<double> energy, energyBeta;
    QList<QString> elements, uElements, symmetry, symmetryBeta;
    QList<QList<double> > mullikenMu, simpleMu, overlapMatrix;
    QList<QList<QList<QList<double> > > > eigenVectors, eigenVectorsBeta, mullikenVectors, mullikenVectorsBeta, lowdinVector, lowdinVectorBeta;

    void SwitchProg(const QString &str, unsigned short &i, unsigned short &a, unsigned short &s);   //Перегруженная функция для поиска в файле
    bool SwitchProg(const QString &str);
    bool FoundItSimple(QTextStream &stream);
    bool FoundItMulliken(QTextStream &stream);
    bool FoundItLowdin(const QString &str, QTextStream &stream);
};

#endif // FOUTPUT_H
