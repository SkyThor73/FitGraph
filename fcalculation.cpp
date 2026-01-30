#include "fcalculation.h"

FCalculation::FCalculation(QObject *parent)
    : QObject{parent}
{}

void FCalculation::Calculation(const FOutput &output, const QList<FElement *> &element)
{
    QList<QList<QList<QList<double> > > > uSVectors = ConvolusionSPD(output, element, output.getOrbital_count(), output.getEigenVectors());   //Создание вектора с вкладами по уникальным атомам
    RefreshVector(output, output.getOrbital_count(), spdSimple);    //Перезаполнение spd вектора
    QList<double> sumSimple = sumSPD(output, output.getOrbital_count(), uSVectors);  //Вектор суммы значений эйгенвекторов по всем атомам
    FillingVector(output, output.getOrbital_count(), uSVectors, sumSimple, spdSimple);  //Заполнение простыми вкладами

    if (output.getUhf()) {    //Выполняет вычисление бета набора если есть uhf
        QList<QList<QList<QList<double> > > > uSVectorsB = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getEigenVectorsBeta());   //Создание вектора с вкладами по уникальным атомам
        RefreshVector(output, output.getOrbital_count_beta(), spdSimpleBeta);    //Перезаполнение spd вектора
        QList<double> sumSimpleBeta = sumSPD(output, output.getOrbital_count_beta(), uSVectorsB);  //Вектор суммы значений эйгенвекторов по всем атомам
        FillingVector(output, output.getOrbital_count_beta(), uSVectorsB, sumSimpleBeta, spdSimpleBeta);  //Заполнение простыми beta вкладами
    }

    if (output.getMulliken()) {    //Выполняет вычисление малликеновских вкладов если NPRINT = 8
        QList<QList<QList<QList<double> > > > uMVectors = ConvolusionSPD(output, element, output.getOrbital_count(), output.getMullikenVectors());   //Создание вектора с вкладами по уникальным атомам
        RefreshVector(output, output.getOrbital_count(), spdMulliken);    //Перезаполнение spd вектора
        QList<double> sumMulliken = sumSPD(output, output.getOrbital_count(), uMVectors);  //Вектор суммы значений векторов по всем атомам
        FillingVector(output, output.getOrbital_count(), uMVectors, sumMulliken, spdMulliken);  //Заполнение млликеновскими вкладами

        if (output.getUhf()) {    //Выполняет вычисление бета набора если есть uhf
            QList<QList<QList<QList<double> > > > uMVectorsB = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getMullikenVectorsBeta());   //Создание вектора с вкладами по уникальным атомам
            RefreshVector(output, output.getOrbital_count_beta(), spdMullikenBeta);    //Перезаполнение spd вектора
            QList<double> sumMullikenBeta = sumSPD(output, output.getOrbital_count_beta(), uMVectorsB);  //Вектор суммы значений векторов по всем атомам
            FillingVector(output, output.getOrbital_count_beta(), uMVectorsB, sumMullikenBeta, spdMullikenBeta);  //Заполнение млликеновскими beta вкладами
        }
    }

    if (output.getLowdin()) {    //Выполняет вычисление вкладов лёвдина если NPRINT = 9
        QList<QList<QList<QList<double> > > > uLVectors = ConvolusionSPD(output, element, output.getOrbital_count(), output.getLowdinVector());   //Создание вектора с вкладами по уникальным атомам
        RefreshVector(output, output.getOrbital_count(), spdLowdin);    //Перезаполнение spd вектора
        QList<double> sumLowdin = sumSPD(output, output.getOrbital_count(), uLVectors);  //Вектор суммы значений векторов по всем атомам
        FillingVector(output, output.getOrbital_count(), uLVectors, sumLowdin, spdLowdin);  //Заполнение лёвдин вкладами

        if (output.getUhf()) {    //Выполняет вычисление бета набора если есть uhf
            QList<QList<QList<QList<double> > > > uLVectorsB = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getLowdinVectorBeta());   //Создание вектора с вкладами по уникальным атомам
            RefreshVector(output, output.getOrbital_count_beta(), spdLowdinBeta);    //Перезаполнение spd вектора
            QList<double> sumLowdinBeta = sumSPD(output, output.getOrbital_count_beta(), uLVectorsB);  //Вектор суммы значений векторов по всем атомам
            FillingVector(output, output.getOrbital_count_beta(), uLVectorsB, sumLowdinBeta, spdLowdinBeta);  //Заполнение лёвдин beta вкладами
        }
    }

    QList<unsigned short> temp_element_number;
    for (int i = 0; i < output.getUnique_atoms(); i++)
        temp_element_number.push_back(element[i]->getElementNumber());
    for (int i = 0; i < output.getUnique_atoms(); i++) {    //Сортировка столбцов пузырьком
        for (int j = 0; j < output.getUnique_atoms() - 1; j++) {
            if (temp_element_number[j] < temp_element_number[j + 1]) {
                temp_element_number.swapItemsAt(j, j + 1);
                for (unsigned short o = 0; o < output.getOrbital_count(); o++)    //Сортировка вектора
                    spdSimple[o].swapItemsAt(j, j + 1);
                if (!spdSimpleBeta.isEmpty()) {
                    for (unsigned short o = 0; o < output.getOrbital_count_beta(); o++)    //Сортировка если есть бета набор
                        spdSimpleBeta[o].swapItemsAt(j, j + 1);
                }
                if (!spdMulliken.isEmpty()) {
                    for (unsigned short o = 0; o < output.getOrbital_count(); o++)    //Сортировка если есть малликен
                        spdMulliken[o].swapItemsAt(j, j + 1);
                    if (!spdMullikenBeta.isEmpty()) {
                        for (unsigned short o = 0; o < output.getOrbital_count_beta(); o++)    //Сортировка если есть бета малликен
                            spdMullikenBeta[o].swapItemsAt(j, j + 1);
                    }
                }
                if (!spdLowdin.isEmpty()) {
                    for (unsigned short o = 0; o < output.getOrbital_count(); o++)    //Сортировка если есть малликен
                        spdLowdin[o].swapItemsAt(j, j + 1);
                    if (!spdLowdinBeta.isEmpty()) {
                        for (unsigned short o = 0; o < output.getOrbital_count_beta(); o++)    //Сортировка если есть бета малликен
                            spdLowdinBeta[o].swapItemsAt(j, j + 1);
                    }
                }
            }
        }
    }
}

void FCalculation::CalculationIntegralArea(const FOutput &output, const FDataPCS &pcs, const unsigned short &set_orbitals, const QList<double> &setEnergy, const QList<QList<QList<QList<double> > > > &spdSet)
{
    //Расчёт высоты отдельных подстолбцов
    //Расчёт интегральной площади
    //Расчёт высоты распределения
    RIntensity.clear();
    RIntensity.resize(set_orbitals);
    integral_area.clear();
    integral_area.resize(set_orbitals);
    betaMO.clear();
    betaMO.resize(set_orbitals);
    RefreshVector(output, set_orbitals, spdRIntensity);
    for (unsigned short o = 0; o < set_orbitals; o++) {
        double betaMO_a = 0.0, betaMO_b = 0.0;
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short spdSetSize = spdSet[o][a][s].size();
                for (unsigned short n = 0; n < spdSetSize; n++) {
                    betaMO_a += spdSet[o][a][s][n] * pcs.getPhotoionizationCrossSection()[source][a][s][n] * pcs.getAsymmetry()[source][a][s][n];
                    betaMO_b += spdSet[o][a][s][n] * pcs.getPhotoionizationCrossSection()[source][a][s][n];
                }
            }
        }
        betaMO[o] = betaMO_a / betaMO_b;
    }
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short spdSetSize = spdSet[o][a][s].size();
                for (unsigned short n = 0; n < spdSetSize; n++) {
                    spdRIntensity[o][a][s][n] += 0.36952117148753 * (2 + betaMO[o] * 0.5) * 2.0 * (spdSet[o][a][s][n] * pcs.getPhotoionizationCrossSection()[source][a][s][n]) / resolution_ratio;   //Первый множитель - высота столбика 0.36952117148753 * 2.0
                    integral_area[o] += ((2 + betaMO[o] * 0.5) * spdSet[o][a][s][n] * pcs.getPhotoionizationCrossSection()[source][a][s][n]);
                    RIntensity[o] += spdRIntensity[o][a][s][n];
                }
            }
        }
    }

    CalculationXY(set_orbitals, setEnergy);
}

void FCalculation::CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy)
{
    abscissaX.clear();
    abscissaX = { initial_value };
    for (unsigned short x = 0; x < steps - 1; x++) {
        abscissaX.push_back(abscissaX[x] + (end_value - initial_value) / steps);
    }

    gaussian.clear();
    gaussian.resize(set_orbitals);
    for (auto& o : gaussian) {
        o.resize(steps);
    }


    for (unsigned short o = 0; o < set_orbitals; o++) {
        if (setEnergy[o] < getUpper_bound()) {
            for (unsigned short x = 0; x < steps; x++) {
                gaussian[o][x] += integral_area[o] * (0.108524845891868 * resolution_ratio / (pow(abscissaX[x] - setEnergy[o], 2) + pow(resolution_ratio, 2) / 4) + 0.304942959407588 / resolution_ratio * pow(M_E, -(2.77258872223978 * pow((abscissaX[x] - setEnergy[o]) / resolution_ratio, 2))));    //Псевдо-Фойгт
            }
        }
    }

    ordinateY.clear();
    for (unsigned short x = 0; x < steps; x++) {
        ordinateY.push_back(0);
        for (unsigned short o = 0; o < set_orbitals; o++) {
            ordinateY[x] += gaussian[o][x];
        }
    }

    if (ordinateY.size() != 0)
        max_value_y = *std::max_element(ordinateY.begin(), ordinateY.end());
}

void FCalculation::CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy, QList<double> resolution_vector)
{
    abscissaX.clear();
    abscissaX = { initial_value };
    for (unsigned short x = 0; x < steps - 1; x++) {
        abscissaX.push_back(abscissaX[x] + (end_value - initial_value) / steps);
    }

    gaussian.clear();
    gaussian.resize(set_orbitals);
    for (auto& o : gaussian) {
        o.resize(steps);
    }

    for (unsigned short o = 0; o < set_orbitals; o++) {
        if (setEnergy[o] < getUpper_bound()) {
            for (unsigned short x = 0; x < steps; x++) {
                gaussian[o][x] += integral_area[o] * (0.108524845891868 * resolution_vector[o] / (pow(abscissaX[x] - setEnergy[o], 2) + pow(resolution_vector[o], 2) / 4) + 0.304942959407588 / resolution_vector[o] * pow(M_E, -(2.77258872223978 * pow((abscissaX[x] - setEnergy[o]) / resolution_vector[o], 2))));    //Псевдо-Фойгт
            }
        }
    }

    ordinateY.clear();
    for (unsigned short x = 0; x < steps; x++) {
        ordinateY.push_back(0);
        for (unsigned short o = 0; o < set_orbitals; o++) {
            ordinateY[x] += gaussian[o][x];
        }
    }

    if (ordinateY.size() != 0)
        max_value_y = *std::max_element(ordinateY.begin(), ordinateY.end());
}

void FCalculation::CalculationEachWidth(const FOutput &output, const FDataPCS &pcs, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &spdSet, QList<double> resolution_vector)
{
    RIntensity.clear();
    RIntensity.resize(set_orbitals);
    RefreshVector(output, set_orbitals, spdRIntensity);
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short spdSetSize = spdSet[o][a][s].size();
                for (unsigned short n = 0; n < spdSetSize; n++) {
                    spdRIntensity[o][a][s][n] += 0.36952117148753 * (2 + betaMO[o] * 0.5) * 2.0 * (spdSet[o][a][s][n] * pcs.getPhotoionizationCrossSection()[source][a][s][n]) / resolution_vector[o];   //Последний множитель - высота столбика
                    RIntensity[o] += spdRIntensity[o][a][s][n];
                }
            }
        }
    }
}

const FCalculation::Sources &FCalculation::getSource() const
{
    return source;
}

void FCalculation::setSource(Sources newSource)
{
    source = newSource;
}

const unsigned short &FCalculation::getSteps() const
{
    return steps;
}

void FCalculation::setSteps(unsigned short newSteps)
{
    steps = newSteps;
}

const double &FCalculation::getInitial_value() const
{
    return initial_value;
}

void FCalculation::setInitial_value(double newInitial_value)
{
    initial_value = newInitial_value;
}

const double &FCalculation::getEnd_value() const
{
    return end_value;
}

void FCalculation::setEnd_value(double newEnd_value)
{
    end_value = newEnd_value;
}

const double &FCalculation::getResolution_ratio() const
{
    return resolution_ratio;
}

void FCalculation::setResolution_ratio(double newResolution_ratio)
{
    resolution_ratio = newResolution_ratio;
}

const double &FCalculation::getMax_value_y() const
{
    return max_value_y;
}

const double &FCalculation::getUpper_bound() const
{
    return upper_bound;
}

void FCalculation::setUpper_bound(double newUpper_bound)
{
    upper_bound = newUpper_bound;
}

const QList<double> &FCalculation::getIntegral_area() const
{
    return integral_area;
}

const QList<double> &FCalculation::getAbscissaX() const
{
    return abscissaX;
}

const QList<double> &FCalculation::getOrdinateY() const
{
    return ordinateY;
}

const QList<double> &FCalculation::getResolution_vector() const
{
    return resolution_vector;
}

void FCalculation::setResolution_vector(const QList<double> &newResolution_vector)
{
    resolution_vector = newResolution_vector;
}

const QList<double> &FCalculation::getRIntensity() const
{
    return RIntensity;
}

const QList<QList<double> > &FCalculation::getGaussian() const
{
    return gaussian;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdSimple() const
{
    return spdSimple;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdSimpleBeta() const
{
    return spdSimpleBeta;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdMulliken() const
{
    return spdMulliken;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdMullikenBeta() const
{
    return spdMullikenBeta;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdLowdin() const
{
    return spdLowdin;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdLowdinBeta() const
{
    return spdLowdinBeta;
}

const QList<QList<QList<QList<double> > > > &FCalculation::getSpdRIntensity() const
{
    return spdRIntensity;
}

void FCalculation::RefreshVector(const FOutput &output, const unsigned short &set_orbitals, QList<QList<QList<QList<double> > > > &spd)
{
    spd.clear();
    spd.resize(set_orbitals);
    for (unsigned short i = 0; i < set_orbitals; i++) {
        spd[i].resize(output.getUnique_atoms());
        for (unsigned short j = 0; j < output.getUnique_atoms(); j++) {
            spd[i][j].resize(output.getKSpdfghi());
            for (unsigned short k = 0; k < output.getKSpdfghi(); k++) {
                spd[i][j][k].resize(19);
            }
        }
    }
}

void FCalculation::FillingVector(const FOutput &output, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &uVectors, const QList<double> &sumVector, QList<QList<QList<QList<double> > > > &spdVectorAB)
{
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short uVectorsSize = uVectors[o][a][s].size();
                for (unsigned short n = 0; n < uVectorsSize; n++) {
                    if (sumVector[o] != 0)
                        spdVectorAB[o][a][s][n] += uVectors[o][a][s][n] / sumVector[o];
                }
            }
        }
    }
}

QList<double> FCalculation::sumSPD(const FOutput &output, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &uVectors)
{
    QList<double> sumVector(set_orbitals);   //Полная сумма значений эйгенвекторов
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short uVectorsSize = uVectors[o][a][s].size();
                for (unsigned short n = 0; n < uVectorsSize; n++)
                    sumVector[o] += uVectors[o][a][s][n];
            }
        }
    }
    return sumVector;
}

QList<QList<QList<QList<double> > > > FCalculation::ConvolusionSPD(const FOutput &output, const QList<FElement *> &element, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &setVectors)
{
    QList<QList<QList<QList<double> > > > sumVectorsSPD(set_orbitals, QList<QList<QList<double> > > (output.getAtom_count(), QList<QList<double> > (output.getKSpdfghi(), QList<double>(19))));    //Суммирование эйгенвекторов по компонентам
    unsigned short i = 0;
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < output.getAtom_count(); a++) {
            i = 0;
            for (const auto &var : output.getUElements()) {
                if (var == output.getElements()[a])
                    break;
                ++i;
            }
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                unsigned short setVectorsSize = setVectors[o][a][s].size();
                switch (s) {
                case 0: {       //S оболочка
                    if (output.getBasis().contains("def2", Qt::CaseInsensitive) || output.getBasis().contains("vp", Qt::CaseInsensitive) || output.getBasis().contains("v-p", Qt::CaseInsensitive)) {
                        if (output.getBasis().contains("s", Qt::CaseInsensitive)) {         //SVP
                            if (element[i]->getElementNumber() < 3) {       //H, He
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 11) {
                                for (unsigned short n = 0; n < 2; n++)      //1s 2s
                                    sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                                for (unsigned short n = 2; n < setVectorsSize; n++)    //в 2s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 19) {
                                for (unsigned short n = 0; n < 3; n++)      //1s 2s 3s
                                    sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                                for (unsigned short n = 3; n < setVectorsSize; n++)    //в 3s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 37) {
                                for (unsigned short n = 0; n < 4; n++)      //1s 2s 3s 4s
                                    sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                                for (unsigned short n = 4; n < setVectorsSize; n++)    //в 4s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 2; n++)      //4s
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                                for (unsigned short n = 2; n < setVectorsSize; n++)    //в 5s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][4] += setVectors[o][a][s][n];
                            }/* else if (element[i]->getElementNumber() < 55) {   //Xe?
                                for (unsigned short n = 0; n < 4; n++)      //4s
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                                for (unsigned short n = 4; n < setVectorsSize; n++)    //в 5s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][4] += setVectors[o][a][s][n];
                            }*/
                        } else if (output.getBasis().contains("t", Qt::CaseInsensitive)) {  //TZVP
                            if (element[i]->getElementNumber() < 3) {       //H, He
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 11) {
                                for (unsigned short n = 0; n < 2; n++)      //1s
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                                for (unsigned short n = 2; n < setVectorsSize; n++)    //в 2s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 19) {
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][0];    //1s
                                for (unsigned short n = 1; n < 3; n++)      //2s
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                                for (unsigned short n = 3; n < setVectorsSize; n++)    //в 3s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 37) {
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][0];    //1s
                                sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][1];    //2s
                                for (unsigned short n = 2; n < 4; n++)      //3s
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                                for (unsigned short n = 4; n < setVectorsSize; n++)    //в 4s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 4; n++)      //4s
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                                for (unsigned short n = 4; n < setVectorsSize; n++)    //в 5s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][4] += setVectors[o][a][s][n];
                            }
                        } else if (output.getBasis().contains("q", Qt::CaseInsensitive)) {  //QZVP
                            if (element[i]->getElementNumber() < 3) {       //H, He
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 11) {   //Ne
                                for (unsigned short n = 0; n < 3; n++)      //1s
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                                for (unsigned short n = 3; n < setVectorsSize; n++)    //в 2s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 19) {   //Ar
                                for (unsigned short n = 0; n < 3; n++)      //1s
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                                for (unsigned short n = 3; n < 6; n++)      //2s
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                                for (unsigned short n = 6; n < setVectorsSize; n++)    //в 3s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 37) {   //Kr
                                for (unsigned short n = 0; n < 3; n++)      //1s
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                                for (unsigned short n = 3; n < 5; n++)      //2s
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                                for (unsigned short n = 5; n < 8; n++)      //3s
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                                for (unsigned short n = 8; n < setVectorsSize; n++)    //в 4s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 4; n++)      //4s
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                                for (unsigned short n = 4; n < setVectorsSize; n++)    //в 5s суммируем оставшиеся коэффициенты
                                    sumVectorsSPD[o][a][s][4] += setVectors[o][a][s][n];
                            }
                        }
                    } else if (output.getBasis().contains("cc", Qt::CaseInsensitive) || output.getBasis().contains("pc", Qt::CaseInsensitive)) {
                        if (element[i]->getElementNumber() < 3) {       //H, He
                            for (unsigned short n = 0; n < setVectorsSize; n++)
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];
                        } else if (element[i]->getElementNumber() < 11) {
                            for (unsigned short n = 0; n < 2; n++)      //1s 2s
                                sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                            for (unsigned short n = 2; n < setVectorsSize; n++)    //в 2s суммируем оставшиеся коэффициенты
                                sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];
                        } else if (element[i]->getElementNumber() < 19) {
                            for (unsigned short n = 0; n < 3; n++)      //1s 2s 3s
                                sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                            for (unsigned short n = 3; n < setVectorsSize; n++)    //в 3s суммируем оставшиеся коэффициенты
                                sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];
                        } else if (element[i]->getElementNumber() < 37) {
                            for (unsigned short n = 0; n < 4; n++)      //1s 2s 3s 4s
                                sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][n];
                            for (unsigned short n = 4; n < setVectorsSize; n++)    //в 4s суммируем оставшиеся коэффициенты
                                sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];
                        }
                    }
                    break;
                }
                case 1: {       //P оболочка
                    if (output.getBasis().contains("def2", Qt::CaseInsensitive) || output.getBasis().contains("vp", Qt::CaseInsensitive) || output.getBasis().contains("v-p", Qt::CaseInsensitive)) {
                        if (output.getBasis().contains("s", Qt::CaseInsensitive)) {             //SVP
                            if (element[i]->getElementNumber() < 3) {           //H, He
                                sumVectorsSPD[o][a][s][0] = 0.0;
                            } else if (element[i]->getElementNumber() < 11) {   //Ne
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                            } else if (element[i]->getElementNumber() < 19) {   //Ar
                                for (unsigned short n = 0; n < 3; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 3; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                            } else if (element[i]->getElementNumber() < 37) {   //Kr
                                for (unsigned short n = 0; n < 3; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 3; n < 6; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                                for (unsigned short n = 6; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 6; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                                for (unsigned short n = 6; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];    //5p
                            }/* else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 9; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                                for (unsigned short n = 9; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];    //5p
                            }*/
                        } else if (output.getBasis().contains("t", Qt::CaseInsensitive)) {      //TZVP
                            if (element[i]->getElementNumber() < 3) {           //H, He
                                sumVectorsSPD[o][a][s][0] = 0.0;
                            } else if (element[i]->getElementNumber() < 11) {   //Ne
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                            } else if (element[i]->getElementNumber() < 19) {   //Ar
                                for (unsigned short n = 0; n < 6; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 6; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                            } else if (element[i]->getElementNumber() < 37) {   //Kr
                                for (unsigned short n = 0; n < 3; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 3; n < 6; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                                for (unsigned short n = 6; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 9; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                                for (unsigned short n = 9; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];    //5p
                            }
                        } else if (output.getBasis().contains("q", Qt::CaseInsensitive)) {      //QZVP
                            if (element[i]->getElementNumber() < 3) {           //H, He
                                sumVectorsSPD[o][a][s][0] = 0.0;
                            } else if (element[i]->getElementNumber() < 11) {   //Ne
                                for (unsigned short n = 0; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                            } else if (element[i]->getElementNumber() < 19) {   //Ar
                                for (unsigned short n = 0; n < 6; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 6; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                            } else if (element[i]->getElementNumber() < 37) {   //Kr
                                for (unsigned short n = 0; n < 3; n++)
                                    sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                                for (unsigned short n = 3; n < 9; n++)
                                    sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                                for (unsigned short n = 9; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                            } else if (element[i]->getElementNumber() < 55) {   //Xe
                                for (unsigned short n = 0; n < 9; n++)
                                    sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                                for (unsigned short n = 9; n < setVectorsSize; n++)
                                    sumVectorsSPD[o][a][s][3] += setVectors[o][a][s][n];    //5p
                            }
                        }
                    } else if (output.getBasis().contains("cc", Qt::CaseInsensitive) || output.getBasis().contains("pc", Qt::CaseInsensitive)) {
                        if (element[i]->getElementNumber() < 3) {           //H, He
                            sumVectorsSPD[o][a][s][0] = 0.0;
                        } else if (element[i]->getElementNumber() < 11) {   //Ne
                            for (unsigned short n = 0; n < setVectorsSize; n++)
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                        } else if (element[i]->getElementNumber() < 19) {   //Ar
                            for (unsigned short n = 0; n < 3; n++)
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                            for (unsigned short n = 3; n < setVectorsSize; n++)
                                sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                        } else if (element[i]->getElementNumber() < 37) {   //Kr
                            for (unsigned short n = 0; n < 3; n++)
                                sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //2p
                            for (unsigned short n = 3; n < 6; n++)
                                sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //3p
                            for (unsigned short n = 6; n < setVectorsSize; n++)
                                sumVectorsSPD[o][a][s][2] += setVectors[o][a][s][n];    //4p
                        }
                    }
                    break;
                }
                case 2: {       //D оболочка
                    if (element[i]->getElementNumber() < 19) {   //Ar
                        sumVectorsSPD[o][a][s][0] = 0.0;
                    } else if (element[i]->getElementNumber() < 37) {   //Kr
                        for (unsigned short n = 0; n < setVectorsSize; n++)
                            sumVectorsSPD[o][a][s][0] += setVectors[o][a][s][n];    //3d
                    } else if (element[i]->getElementNumber() < 55) {   //Xe
                        for (unsigned short n = 0; n < setVectorsSize; n++)
                            sumVectorsSPD[o][a][s][1] += setVectors[o][a][s][n];    //4d
                    }
                    break;
                }
                case 3: {       //F оболочка
                    if (output.getProgramm() == FOutput::GAUSSIAN || output.getProgramm() == FOutput::ORCA) {
                        for (unsigned short n = 0; n < setVectorsSize / 7; n++) {
                            for (unsigned short i = (n * 7); i < ((n * 7) + 7); i++) {
                                sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][i];
                            }
                        }
                    } else {
                        for (unsigned short n = 0; n < setVectorsSize / 10; n++) {
                            for (unsigned short i = (n * 10); i < ((n * 10) + 10); i++) {
                                sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][i];
                            }
                        }
                    }
                    break;
                }
                case 4: {       //G оболочка
                    for (unsigned short n = 0; n < setVectorsSize / 15; n++) {
                        for (unsigned short i = (n * 15); i < ((n * 15) + 15); i++) {
                            sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][i];
                        }
                    }
                    break;
                }
                case 5: {       //H оболочка
                    for (unsigned short n = 0; n < setVectorsSize / 21; n++) {
                        for (unsigned short i = (n * 21); i < ((n * 21) + 21); i++) {
                            sumVectorsSPD[o][a][s][n] += setVectors[o][a][s][i];
                        }
                    }
                    break;
                }
                default:
                    break;
                }
            }
        }
    }

    QList<QList<QList<QList<double> > > > uVectors(set_orbitals, QList<QList<QList<double> > > (output.getUnique_atoms(), QList<QList<double> > (output.getKSpdfghi(), QList<double>(19))));
    for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
        for (unsigned short aa = 0; aa < output.getAtom_count(); aa++) {
            if (output.getUElements()[a] == output.getElements()[aa]) { //Суммирование эйгенвекторов по одинаковым атомам
                for (unsigned short o = 0; o < set_orbitals; o++) {
                    for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                        unsigned short sumVectorsSPDSize = sumVectorsSPD[o][a][s].size();
                        for (unsigned short n = 0; n < sumVectorsSPDSize; n++)
                            uVectors[o][a][s][n] += sumVectorsSPD[o][aa][s][n];
                    }
                }
            }
        }
    }
    return uVectors;
}
