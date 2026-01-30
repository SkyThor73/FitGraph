#include "fadjusting.h"

FAdjusting::FAdjusting(QObject *parent)
    : QObject{parent}
{}

void FAdjusting::FindBetterShift(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp)
{
    unsigned short shift_steps = 0;
    double step_energy = (calculate.getEnd_value() - calculate.getInitial_value()) / calculate.getSteps(); //Энергетический шаг
    switch (calculate.getSource()) {
    case FCalculation::HE_I:
        shift_steps = 3.0  / step_energy;       //Количество точек на которое необходимо сдвинуть расчёт относительно эксперимента
        if (shift_steps > calculate.getSteps())
            shift_steps = calculate.getSteps() - 1;
        break;
    case FCalculation::HE_II:
        shift_steps = 5.0  / step_energy;       //Количество точек на которое необходимо сдвинуть расчёт относительно эксперимента
        if (shift_steps > calculate.getSteps())
            shift_steps = calculate.getSteps() - 1;
        break;
    case FCalculation::MG:
        shift_steps = 5.0  / step_energy;       //Количество точек на которое необходимо сдвинуть расчёт относительно эксперимента
        fwhm_max = 5;
        if (shift_steps > calculate.getSteps())
            shift_steps = calculate.getSteps() - 1;
        break;
    case FCalculation::AL:
        shift_steps = 5.0  / step_energy;      //Количество точек на которое необходимо сдвинуть расчёт относительно эксперимента
        fwhm_max = 5;
        if (shift_steps > calculate.getSteps())
            shift_steps = calculate.getSteps() - 1;
        break;
    default:
        break;
    }
    //Расширение области построения
    calculate.setInitial_value(calculate.getInitial_value() - shift_steps * step_energy);
    calculate.setEnd_value(calculate.getEnd_value() + shift_steps * step_energy);
    calculate.setSteps(calculate.getSteps() + 2 * shift_steps);
    calculate.CalculationXY(set_orbitals, setEnergy);
    QList<double> calculatedY = calculate.getOrdinateY(); //Временный набор расчётных игреков
    QList<double> experimentalY = exp.getExpOrdinateY();  //Временный набор экспериментальных игреков
    double multiplier = calculate.getMax_value_y() / exp.getMax_value_exp_y();
    for (int i = 0; i < exp.getExpOrdinateY().size(); ++i)
        experimentalY[i] *= multiplier;
    QList<double> squareInterval;
    double sum_value_interval = 0.0;
    for(unsigned short i = 0; i != (2 * shift_steps); ++i) {   //Сдвиг экспериментального относительно расчётного по всему расчётному диапазону
        sum_value_interval = 0.0;
        for(unsigned short j = 0; j < experimentalY.size(); ++j)
            sum_value_interval += pow(calculatedY[i + j] - experimentalY[j], 2.0);  //МНК
        for (unsigned short k = 0; k < i; ++k)
            sum_value_interval += pow(calculatedY[k], 2.0);
        for (unsigned short k = (i + experimentalY.size()); k < calculatedY.size(); ++k)
            sum_value_interval += pow(calculatedY[k], 2.0);
        squareInterval.push_back(sum_value_interval);   //Заполнение вектора значениями площадей
    }
    unsigned short iterator = std::distance(squareInterval.begin(), std::min_element(squareInterval.begin(), squareInterval.end()));   //Итератор в котором содержится минимальное значение
    shift_x = step_energy * (shift_steps - iterator);   //Значение на которое сдвинутся расчётные уровни
    visible_shift_x += shift_x; //Отображаемое значение сдвига
    //Возврат области построения к исходной
    calculate.setInitial_value(calculate.getInitial_value() + shift_steps * step_energy);
    calculate.setEnd_value(calculate.getEnd_value() - shift_steps * step_energy);
    calculate.setSteps(calculate.getSteps() - 2 * shift_steps);
}

void FAdjusting::FindBetterVise(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp)
{
    // QList<double> temp_resolution_vector = calculate.getResolution_vector();
    // if(temp_resolution_vector.size() == 5) {
    //     temp_resolution_vector[0] = 1.036;
    //     temp_resolution_vector[1] = 3.26;
    //     temp_resolution_vector[2] = 1.27;
    //     temp_resolution_vector[3] = 1.67;
    //     temp_resolution_vector[4] = 0.88;
    // }
    // calculate.CalculationXY(set_orbitals, setEnergy, temp_resolution_vector);
    // calculate.setResolution_vector(temp_resolution_vector);
    double temp_sum_exp = 0.0, temp_sum_cal = 0.0;
    for (auto var : exp.getExpOrdinateY())
        temp_sum_exp += var;
    for (auto var : calculate.getOrdinateY())
        temp_sum_cal += var;
    double multiplier = temp_sum_cal / temp_sum_exp;        //Множитель точек нормированный к площади
    //Приравнивание интегральной площади
    QList<double> experimentalY = exp.getExpOrdinateY();
    for (int i = 0; i < exp.getExpOrdinateY().size(); ++i)
        experimentalY[i] *= multiplier;
    exp.setExpOrdinateY(experimentalY);
    exp.setMax_value_exp_y(*std::max_element(experimentalY.begin(), experimentalY.end()));  //Максимальная высота может измениться
}

void FAdjusting::FindBetterWidth(const unsigned short &set_orbitals, const QList<double> &setEnergy, FCalculation &calculate, FExperimental &exp)
{
    QList<double> squareInterval;     //Вектор содержащий в себе средние значения интервалов
    squareInterval.push_back(specific_min_square);  //Последнее значение соответствия
    double decimal = 0.0;               //Порядок точности слагаемого
    unsigned short same_value = 0;              //Определение количества орбиталей с одинаковой энергией подряд

    QList<double> temp_resolution_vector = calculate.getResolution_vector();

    for (unsigned short o = 0; o < temp_resolution_vector.size(); ++o) {    //Начало с последней орбитали
        if (setEnergy[o] > (calculate.getInitial_value() - 2) && setEnergy[o] < (calculate.getEnd_value() + 2) && setEnergy[o] < calculate.getUpper_bound()) {
            same_value = 0;              //Определение количества орбиталей с одинаковой энергией подряд
            while ((o + same_value + 1) != setEnergy.size() && setEnergy[o] == setEnergy[o + same_value + 1]) {
                ++same_value;
            }
            decimal = 0.9;
            squareInterval.clear();
            squareInterval.push_back(specific_min_square);  //Последнее значение соответствия
            for (unsigned short j = 1; ; j++) {
                for (unsigned short s = 0; s < same_value + 1; ++s) {   //Изменение ширины полос всех повторяющихся МО
                    temp_resolution_vector[o + s] += (decimal / (same_value + 1));
                }
                if (temp_resolution_vector[o] >= fwhm_min && temp_resolution_vector[o] <= fwhm_max) {
                    calculate.CalculationXY(set_orbitals, setEnergy, temp_resolution_vector);
                    CalculationCompliance(calculate, exp);
                    squareInterval.push_back(specific_min_square);
                    if (squareInterval[j] < squareInterval[j - 1]) {
                        if ((decimal > 0 && decimal < 0.0001) || (decimal > -0.0001 && decimal < 0)) {
                            break;
                        }
                        decimal *= -0.1;
                    }
                } else {
                    for (unsigned short s = 0; s < same_value + 1; ++s) {   //Изменение ширины полос всех повторяющихся МО
                        temp_resolution_vector[o + s] -= (decimal / (same_value + 1));
                    }
                    break;
                }
            }
            o += same_value;
        }
    }
    calculate.setResolution_vector(temp_resolution_vector);
}

void FAdjusting::ShiftEnergy(QList<double> &new_energy)
{
    for (unsigned short o = 0; o < new_energy.size(); o++)
        new_energy[o] += shift_x;
}

void FAdjusting::CalculationCompliance(const FCalculation &calculate, const FExperimental &exp)
{
    double sum_exp = 0.0, sum_square_res = 0.0, sum_square_tot = 0.0, mean_exp = 0.0;
    for (auto var : exp.getExpOrdinateY())
        sum_exp += var;
    mean_exp = sum_exp / calculate.getSteps();
    for (int i = 0; i < exp.getExpOrdinateY().size(); ++i)
        sum_square_res += pow(exp.getExpOrdinateY()[i] - calculate.getOrdinateY()[i], 2.0);
    for (auto var : exp.getExpOrdinateY())
        sum_square_tot += pow((var - mean_exp), 2.0);
    specific_min_square = 1 - sum_square_res / sum_square_tot;
}

const double &FAdjusting::getSpecific_min_square() const
{
    return specific_min_square;
}

const double &FAdjusting::getVisible_shift_x() const
{
    return visible_shift_x;
}

void FAdjusting::setVisible_shift_x(double newVisible_shift_x)
{
    visible_shift_x = newVisible_shift_x;
}
