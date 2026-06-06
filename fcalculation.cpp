#include "fcalculation.h"

FCalculation::FCalculation(QObject *parent)
    : QObject{parent}
{}

void FCalculation::Calculation(const FOutput &output, const QList<FElement *> &element)
{
    //Базовые вклады (Simple)
    spdSimple = ConvolusionSPD(output, element, output.getOrbital_count(), output.getEigenVectors());
    if (output.getUhf()) {
        spdSimpleBeta = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getEigenVectorsBeta());
    }
    //Вклады по Малликену
    if (output.getMulliken()) {
        spdMulliken = ConvolusionSPD(output, element, output.getOrbital_count(), output.getMullikenVectors());
        if (output.getUhf()) {
            spdMullikenBeta = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getMullikenVectorsBeta());
        }
    }
    //Вклады по Лёвдину
    if (output.getLowdin()) {
        spdLowdin = ConvolusionSPD(output, element, output.getOrbital_count(), output.getLowdinVector());
        if (output.getUhf()) {
            spdLowdinBeta = ConvolusionSPD(output, element, output.getOrbital_count_beta(), output.getLowdinVectorBeta());
        }
    }
}

void FCalculation::CalculationIntegralArea(const FOutput &output, const FDataPCS &pcs, const unsigned short &set_orbitals, const QList<double> &setEnergy, const QList<QList<QList<QList<double> > > > &spdSet)
{
    //Расчёт высоты отдельных подстолбцов
    //Расчёт интегральной площади
    //Расчёт высоты распределения
    RIntensity.assign(set_orbitals, 0.0);
    integral_area.assign(set_orbitals, 0.0);
    betaMO.assign(set_orbitals, 0.0);
    //Выделение памяти для spdRIntensity
    spdRIntensity.clear();
    spdRIntensity.resize(set_orbitals);
    for (unsigned short o = 0; o < set_orbitals; ++o) {
        spdRIntensity[o].resize(output.getUnique_atoms());
        for (unsigned short a = 0; a < output.getUnique_atoms(); ++a) {
            spdRIntensity[o][a].resize(output.getKSpdfghi());
            for (unsigned short s = 0; s < output.getKSpdfghi(); ++s) {
                spdRIntensity[o][a][s].assign(spdSet[o][a][s].size(), 0.0);
            }
        }
    }
    for (unsigned short o = 0; o < set_orbitals; o++) {
        double betaMO_a = 0.0, betaMO_b = 0.0;
        // Первый проход: вычисление коэффициента асимметрии betaMO
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                const auto& shell = spdSet[o][a][s];
                for (unsigned short n = 0; n < shell.size(); n++) {
                    const double pop = shell[n];
                    if (pop == 0.0) continue; //Пропуск пустых математических узлов
                    const double cross_section = pcs.getPhotoionizationCrossSection()[source][a][s][n];
                    const double asymmetry = pcs.getAsymmetry()[source][a][s][n];
                    betaMO_a += pop * cross_section * asymmetry;
                    betaMO_b += pop * cross_section;
                }
            }
        }
        //Защита от деления на ноль (для пустых виртуальных орбиталей)
        betaMO[o] = (betaMO_b != 0.0) ? (betaMO_a / betaMO_b) : 0.0;
        const double beta_factor = 2.0 + betaMO[o] * 0.5;
        // Второй проход: расчет интегральной площади
        for (unsigned short a = 0; a < output.getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output.getKSpdfghi(); s++) {
                const auto& shell = spdSet[o][a][s];
                for (unsigned short n = 0; n < shell.size(); n++) {
                    const double pop = shell[n];
                    if (pop == 0.0) continue;
                    const double cross_section = pcs.getPhotoionizationCrossSection()[source][a][s][n];
                    const double intensity_contrib = pop * cross_section;
                    spdRIntensity[o][a][s][n] = 0.36952117148753 * beta_factor * 2.0 * intensity_contrib / resolution_ratio;
                    integral_area[o] += beta_factor * intensity_contrib;
                    RIntensity[o] += spdRIntensity[o][a][s][n];
                }
            }
        }
    }
    CalculationXY(set_orbitals, setEnergy);
}

void FCalculation::CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy)
{
    if (steps == 0) return;
    //Формирование сетки оси X без накопления ошибки округления
    abscissaX.resize(steps);
    const double step_size = (end_value - initial_value) / steps;
    for (unsigned short x = 0; x < steps; x++) {
        abscissaX[x] = initial_value + x * step_size;
    }
    //Аллокация и очистка массивов
    gaussian.resize(set_orbitals);
    for (auto& o : gaussian) {
        o.assign(steps, 0.0);
    }
    ordinateY.assign(steps, 0.0);
    //Предварительный расчет инвариантных констант профиля Псевдо-Фойгта
    const double res2 = resolution_ratio * resolution_ratio;
    const double L_coeff = 0.108524845891868 * resolution_ratio;
    const double G_coeff = 0.304942959407588 / resolution_ratio;
    const double exp_coeff = -2.77258872223978 / res2;
    for (unsigned short o = 0; o < set_orbitals; o++) {
        //Пропускаем расчет, если энергия вне диапазона или площадь равна нулю
        if (setEnergy[o] >= getUpper_bound() || integral_area[o] == 0.0) continue;
        const double E_o = setEnergy[o];
        const double area = integral_area[o];
        for (unsigned short x = 0; x < steps; x++) {
            const double dE = abscissaX[x] - E_o;
            const double dE2 = dE * dE;
            // Расчет компоненты Лоренца
            const double L_val = L_coeff / (dE2 + res2 * 0.25);
            // Расчет компоненты Гаусса
            const double G_val = G_coeff * std::exp(dE2 * exp_coeff);
            const double val = area * (L_val + G_val);
            gaussian[o][x] = val;
            //Прямая аккумуляция в итоговый спектр
            ordinateY[x] += val;
        }
    }
    if (!ordinateY.isEmpty())
        max_value_y = *std::max_element(ordinateY.begin(), ordinateY.end());
}

void FCalculation::CalculationXY(const unsigned short &set_orbitals, const QList<double> &setEnergy, QList<double> resolution_vector)
{
    if (steps == 0) return;
    abscissaX.resize(steps);
    const double step_size = (end_value - initial_value) / steps;
    for (unsigned short x = 0; x < steps; x++) {
        abscissaX[x] = initial_value + x * step_size;
    }
    gaussian.resize(set_orbitals);
    for (auto& o : gaussian) {
        o.assign(steps, 0.0);
    }
    ordinateY.assign(steps, 0.0);
    for (unsigned short o = 0; o < set_orbitals; o++) {
        if (setEnergy[o] >= getUpper_bound() || integral_area[o] == 0.0) continue;
        const double E_o = setEnergy[o];
        const double area = integral_area[o];
        // Инварианты переносятся сюда, так как resolution индивидуально для орбитали
        const double res = resolution_vector[o];
        const double res2 = res * res;
        const double L_coeff = 0.108524845891868 * res;
        const double G_coeff = 0.304942959407588 / res;
        const double exp_coeff = -2.77258872223978 / res2;
        for (unsigned short x = 0; x < steps; x++) {
            const double dE = abscissaX[x] - E_o;
            const double dE2 = dE * dE;
            const double L_val = L_coeff / (dE2 + res2 * 0.25);
            const double G_val = G_coeff * std::exp(dE2 * exp_coeff);
            const double val = area * (L_val + G_val);
            gaussian[o][x] = val;
            ordinateY[x] += val;
        }
    }
    if (!ordinateY.isEmpty())
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

QList<QList<QList<QList<double> > > > FCalculation::ConvolusionSPD(const FOutput &output, const QList<FElement *> &element, const unsigned short &set_orbitals, const QList<QList<QList<QList<double> > > > &setVectors)
{
    //Кэширование типа базисного набора
    enum class BasisFamily { SVP, TZVP, QZVP, DUNNING, UNKNOWN }; //Семейство базисов
    BasisFamily currentBasis = BasisFamily::UNKNOWN;
    const QString basisStr = output.getBasis();
    if (basisStr.contains(u"def2", Qt::CaseInsensitive) || basisStr.contains(u"vp", Qt::CaseInsensitive) || basisStr.contains(u"v-p", Qt::CaseInsensitive)) {
        if (basisStr.contains(u"s", Qt::CaseInsensitive)) currentBasis = BasisFamily::SVP;
        else if (basisStr.contains(u"t", Qt::CaseInsensitive)) currentBasis = BasisFamily::TZVP;
        else if (basisStr.contains(u"q", Qt::CaseInsensitive)) currentBasis = BasisFamily::QZVP;
    } else if (basisStr.contains(u"cc", Qt::CaseInsensitive) || basisStr.contains(u"pc", Qt::CaseInsensitive)) {
        currentBasis = BasisFamily::DUNNING;
    }
    const unsigned short atom_count = output.getAtom_count();
    const unsigned short unique_count = output.getUnique_atoms();
    const unsigned short subshells = output.getKSpdfghi();
    //Хэширование Z-номеров для доступа за O(1)
    QHash<QString, unsigned short> nameToZ;
    for (unsigned short i = 0; i < unique_count; ++i)
        nameToZ.insert(output.getUElements()[i], element[i]->getElementNumber());
    //Создание и сортировка массива строк c названиями атомов в порядке уменьшения зарядового числа
    QList<QString> sorted_elements = output.getUElements();
    std::sort(sorted_elements.begin(), sorted_elements.end(), [&](const QString& a, const QString& b) {
        return nameToZ[a] > nameToZ[b];
    });
    std::vector<unsigned short> unique_map(atom_count, 0);
    std::vector<unsigned short> Z_numbers(unique_count, 0);
    for (unsigned short a = 0; a < atom_count; ++a) {
        for (unsigned short u = 0; u < unique_count; ++u) {
            if (output.getElements()[a] == sorted_elements[u]) {
                unique_map[a] = u;
                Z_numbers[u] = nameToZ[sorted_elements[u]]; // Извлечение Z-номера из хэш-таблицы
                break;
            }
        }
    }
    //Вектор вкладов
    QList<QList<QList<QList<double>>>> uVectors(set_orbitals, QList<QList<QList<double>>>(unique_count, QList<QList<double>>(subshells, QList<double>(19, 0.0))));
    //Свертка неравномерных подоболочек (s, p, d) через границы индексов
    auto contract_shell = [](const auto& src, auto& tgt, auto start_idx, std::initializer_list<unsigned short> bounds) {
        auto tgt_idx = start_idx;
        auto it = bounds.begin();
        for (auto i = 0; i < src.size(); ++i) {
            while (it != bounds.end() && i >= *it) {
                tgt_idx++;
                it++;
            }
            if (tgt_idx < tgt.size()) tgt[tgt_idx] += src[i];
        }
    };
    //Равномерная аккумуляция гармоник (f, g, h)
    auto reduce_angular = [](const auto& src, auto& tgt, auto chunk_size) {
        if (chunk_size == 0) return;
        for (auto i = 0; i < src.size(); ++i) {
            auto tgt_idx = i / chunk_size;
            if (tgt_idx < tgt.size()) tgt[tgt_idx] += src[i];
        }
    };
    for (unsigned short o = 0; o < set_orbitals; o++) {
        for (unsigned short a = 0; a < atom_count; a++) {
            const unsigned short u = unique_map[a];
            const unsigned short Z = Z_numbers[u];
            for (unsigned short s = 0; s < subshells; s++) {
                const auto& source = setVectors[o][a][s];
                auto& target = uVectors[o][u][s];
                if (source.isEmpty()) continue;
                switch (s) {
                case 0: // S-оболочка
                    if (currentBasis == BasisFamily::SVP) {
                        if (Z < 3)       contract_shell(source, target, 0, {});
                        else if (Z < 11) contract_shell(source, target, 0, {1});
                        else if (Z < 19) contract_shell(source, target, 0, {1, 2});
                        else if (Z < 37) contract_shell(source, target, 0, {1, 2, 3});
                        else if (Z < 55) contract_shell(source, target, 3, {2});
                    } else if (currentBasis == BasisFamily::TZVP) {
                        if (Z < 3)       contract_shell(source, target, 0, {});
                        else if (Z < 11) contract_shell(source, target, 0, {2});
                        else if (Z < 19) contract_shell(source, target, 0, {1, 3});
                        else if (Z < 37) contract_shell(source, target, 0, {1, 2, 4});
                        else if (Z < 55) contract_shell(source, target, 3, {4});
                    } else if (currentBasis == BasisFamily::QZVP) {
                        if (Z < 3)       contract_shell(source, target, 0, {});
                        else if (Z < 11) contract_shell(source, target, 0, {3});
                        else if (Z < 19) contract_shell(source, target, 0, {3, 6});
                        else if (Z < 37) contract_shell(source, target, 0, {3, 5, 8});
                        else if (Z < 55) contract_shell(source, target, 3, {4});
                    } else if (currentBasis == BasisFamily::DUNNING) {
                        if (Z < 3)       contract_shell(source, target, 0, {});
                        else if (Z < 11) contract_shell(source, target, 0, {2});
                        else if (Z < 19) contract_shell(source, target, 0, {3});
                        else if (Z < 37) contract_shell(source, target, 0, {4});
                    }
                    break;
                case 1: // P-оболочка
                    if (Z < 3) break; // H, He игнорируются
                    if (currentBasis == BasisFamily::SVP) {
                        if (Z < 11)      contract_shell(source, target, 0, {});
                        else if (Z < 19) contract_shell(source, target, 0, {3});
                        else if (Z < 37) contract_shell(source, target, 0, {3, 6});
                        else if (Z < 55) contract_shell(source, target, 2, {6});
                    } else if (currentBasis == BasisFamily::TZVP) {
                        if (Z < 11)      contract_shell(source, target, 0, {});
                        else if (Z < 19) contract_shell(source, target, 0, {6});
                        else if (Z < 37) contract_shell(source, target, 0, {3, 6});
                        else if (Z < 55) contract_shell(source, target, 2, {9});
                    } else if (currentBasis == BasisFamily::QZVP) {
                        if (Z < 11)      contract_shell(source, target, 0, {});
                        else if (Z < 19) contract_shell(source, target, 0, {6});
                        else if (Z < 37) contract_shell(source, target, 0, {3, 9});
                        else if (Z < 55) contract_shell(source, target, 2, {9});
                    } else if (currentBasis == BasisFamily::DUNNING) {
                        if (Z < 11)      contract_shell(source, target, 0, {});
                        else if (Z < 19) contract_shell(source, target, 0, {3});
                        else if (Z < 37) contract_shell(source, target, 0, {3, 6});
                    }
                    break;
                case 2: // D-оболочка
                    if (Z < 19) break;
                    if (Z < 37)      contract_shell(source, target, 0, {});
                    else if (Z < 55) contract_shell(source, target, 1, {});
                    break;
                case 3: // F-оболочка
                    reduce_angular(source, target, (output.getProgramm() == FOutput::GAUSSIAN || output.getProgramm() == FOutput::ORCA) ? 7 : 10);
                    break;
                case 4: // G-оболочка
                    reduce_angular(source, target, 15);
                    break;
                case 5: // H-оболочка
                    reduce_angular(source, target, 21);
                    break;
                }
            }
        }
    }
    //Финальная физическая редукция заселенностей
    for (auto& orbital : uVectors) {
        for (auto& atom : orbital) {
            for (auto& shell : atom) {
                if (shell.isEmpty()) continue; //Пропуск пустых данных
                double sum = 0.0;
                double max_abs_val = -1.0;
                unsigned short dominant_idx = 0;
                bool has_population = false;
                //Определение доминирующего физического подуровня и алгебраической суммы
                for (auto n = 0; n < shell.size(); ++n) {
                    const double val = shell[n];
                    const double abs_val = std::abs(val);
                    sum += val;
                    if (abs_val > max_abs_val) {
                        max_abs_val = abs_val;
                        dominant_idx = n;
                        if (abs_val > 1e-12) has_population = true; //Фильтр вычислительного шума
                    }
                }
                //Концентрация электронной плотности в нужный подуровень
                if (has_population) {
                    std::fill(shell.begin(), shell.end(), 0.0);
                    shell[dominant_idx] = sum;
                }
            }
        }
    }
    return uVectors;
}
