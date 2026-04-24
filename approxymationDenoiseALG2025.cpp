#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iostream>
#include <string> 
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace chrono;

// Параметры NLM
const int PATCH_RADIUS = 3;    // Размер патча (7x7)
const int SEARCH_RADIUS = 20;   // Область поиска (41x41)
const float H = 10.0f;          // Параметр сглаживания

// Параметр ApproxiP2EPEN
const int THRESHOLD_VALUE_LAPLACIAN = 1000;

//Вычисление дисперсии окна
float Dispersion(vector<unsigned char>& window) {
    float sum_square_deviation = 0.0, median_value, window_sum = 0.0;
    for (int k = 0; k < 25; k++) {
        window_sum += window[k];
    }
    median_value = window_sum / 25;
    for (int k = 0; k < 25; k++) {
        sum_square_deviation += pow(window[k] - median_value, 2);
    }

    return sum_square_deviation / 25;
}

//Вычисление магнитуды градиента окна
float GradientMagnitude(vector<unsigned char>& window) {
    //Маски оператора Собеля
    vector<int> Gx = { -1, -2, 0, 2, 1, -2, -4, 0, 4, 2, -3, -6, 0, 6, 3, -2, -4, 0, 4, 2, -1, -2, 0, 2, 1 };   //координата X
    vector<int> Gy = { -1, -2, -3, -2, -1, -2, -4, -6, -4, -2, 0, 0, 0, 0, 0, 2, 4, 6, 4, 2, 1, 2, 3, 2, 1 };   //координата Y

    float df_gx = 0, df_gy = 0;
    for (int k = 0; k < 25; k++) {
        df_gx += window[k] * Gx[k];
        df_gy += window[k] * Gy[k];
    }

    return sqrt(pow(df_gx, 2) + pow(df_gy, 2));
}

//Вычисление лапласиана окна
float Laplacian(vector<unsigned char>& window) {
    //Маска оператора Лапласа для окна 5 на 5
    vector<int> L = { 0, 0, -1, 0, 0, 0, -1, -2, -1, 0, -1, -2, 16, -2, -1, 0, -1, -2, -1, 0, 0, 0, -1, 0, 0 };
    float laplacian = 0.0;
    for (int k = 0; k < 25; k++) {
        laplacian += window[k] * L[k];
    }
    //Возвращаем модуль лапласиана
    return abs(laplacian);
}

//Функция вычисления центрального значения путём аппроксимации по полиному 2-й степени
int Polynomial2Approximation(vector<unsigned char>& window) {
    //Массивы с возможными значениями координат
    vector<int> xCoords = { -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0, 1, 2 };
    vector<int> yCoords = { 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -2, -2, -2, -2, -2 };
    int sum_values = 0, sum_mult_valueX2 = 0, sum_mult_valueY2 = 0;
    for (int k = 0; k < 25; k++) {
        sum_values += window[k];
        sum_mult_valueX2 += window[k] * pow(xCoords[k], 2);
        sum_mult_valueY2 += window[k] * pow(yCoords[k], 2);
    }
    //Вычисляем значение функции в точке (0, 0)
    float new_central_value = (27 * sum_values - 5 * (sum_mult_valueX2 + sum_mult_valueY2)) / 175;
    int z_ncv = new_central_value;
    //Проверка на попадаие в допустимый интервал значений
    if (new_central_value < 0 || new_central_value > 255) {
        return (new_central_value < 0) ? 0 : 255;
    }
    //Возвращаем округлённое до целого значение
    return (new_central_value - z_ncv < 0.5) ? z_ncv : z_ncv + 1;
}

//Функция вычисления центрального значения путём аппроксимации по положительно-направленной экспоненте
int PositiveExponentialApproximation(vector<unsigned char>& window) {
    float logaryfmic_sum = 0.0;
    for (int k = 0; k < 25; k++) {
        logaryfmic_sum += log(window[k] + 1);
    }
    //Вычисляем значение функции в точке (0, 0)
    float new_central_value = exp(logaryfmic_sum / 25) - 1;
    int z_ncv = new_central_value;

    return (new_central_value - z_ncv < 0.5) ? z_ncv : z_ncv + 1;
}

//Функция вычисления центрального значения путём аппроксимации по отрицательно-направленной экспоненте
int NegotiveExponentialApproximation(vector<unsigned char>& window) { 
    float logaryfmic_sum = 0.0;
    for (int k = 0; k < 25; k++) {
        logaryfmic_sum += log(256 - window[k]);
    }
    //Вычисляем значение функции в точке (0, 0)
    float new_central_value = 256 - exp(logaryfmic_sum / 25);
    int z_ncv = new_central_value;

    return (new_central_value - z_ncv < 0.5) ? z_ncv : z_ncv + 1;
}

//Вычисление лапласиана изменённого окна
float LaplacianNew(vector<unsigned char>& window, int central_value) {
    vector<unsigned char> new_window = window;
    new_window[12] = central_value;     //Изменяем центральный пиксель окна
    //Маска оператора Лапласа для окна 5 на 5
    vector<int> L = { 0, 0, -1, 0, 0, 0, -1, -2, -1, 0, -1, -2, 16, -2, -1, 0, -1, -2, -1, 0, 0, 0, -1, 0, 0 };
    float laplacian = 0;
    for (int k = 0; k < 25; k++) {
        laplacian += new_window[k] * L[k];
    }
    //Возвращаем модуль лапласиана
    return abs(laplacian);
}

//Фенкция для получения нового значения центрального пикселя окна
int ChangeCentralPixel(vector<unsigned char>& window) {
    float dispersion = Dispersion(window);
    float gradient_magnitude = GradientMagnitude(window);
    float laplasian = Laplacian(window);

    //cout << "Дисперсия: " << dispersion << endl;
    //cout << "Магнитуда градиента: " << gradient_magnitude << endl;
    //cout << "Лапласиан: " << laplasian << endl;

    //Возвращаем старое значение пикселя, если шума нет
    if (laplasian < THRESHOLD_VALUE_LAPLACIAN) {
        return window[12];
    }
    //При гладкой поверхности
    if (dispersion <= 256) {
        return Polynomial2Approximation(window);
    }
    //На резких текстурных переходах
    else if (gradient_magnitude > 4599) {
        int new_value_pos_exp = PositiveExponentialApproximation(window);
        int new_value_neg_exp = NegotiveExponentialApproximation(window);
        float laplacian_pos_exp = LaplacianNew(window, new_value_pos_exp);
        float laplacian_neg_exp = LaplacianNew(window, new_value_neg_exp);
        return (laplacian_pos_exp < laplacian_neg_exp) ? new_value_pos_exp : new_value_neg_exp;
    }
    //В остальных случаях (шум, слабые переходы, текстуры)
    else {
        int new_value_polynom = Polynomial2Approximation(window);
        int new_value_pos_exp = PositiveExponentialApproximation(window);
        int new_value_neg_exp = NegotiveExponentialApproximation(window);
        float laplacian_polynom = LaplacianNew(window, new_value_polynom);
        float laplacian_pos_exp = LaplacianNew(window, new_value_pos_exp);
        float laplacian_neg_exp = LaplacianNew(window, new_value_neg_exp);
        if (laplacian_polynom < laplacian_pos_exp && laplacian_polynom < laplacian_neg_exp) {
            return new_value_polynom;
        } else {
            return (laplacian_pos_exp < laplacian_neg_exp) ? new_value_pos_exp : new_value_neg_exp;
        }
    }
    return 255; // По умолчанию
}

//Нахождение медианы окна
int MedianPixel(vector<unsigned char>& window) { 
    sort(begin(window), end(window));
    return window[4];
}

//Нахождение среднего значения окна
int AveragePixel(vector<unsigned char>& window) {
    int sum_pixel_values = 0;
    for (int k = 0; k < 9; k++) {
        sum_pixel_values += window[k];
    }
    return sum_pixel_values / 9;
}

//Нахождение нового значения пикселя окна путём гауссовской фильтрации
int GaussPixel(vector<unsigned char>& window, vector<double>& kernel) {
    double gauss_value = 0.0;
    for (int k = 0; k < window.size(); k++) {
        gauss_value += window[k] * kernel[k];
    }
    int z_gv = gauss_value;

    return (gauss_value - z_gv < 0.5) ? z_gv : z_gv + 1;
}

// Вычисление расстояния между патчами для NLM
float PatchDistance(const unsigned char* img, int width, int height, int channels, int x1, int y1, int x2, int y2) {
    float sum = 0.0f;
    for (int dy = -PATCH_RADIUS; dy <= PATCH_RADIUS; ++dy) {
        for (int dx = -PATCH_RADIUS; dx <= PATCH_RADIUS; ++dx) {
            int px1 = x1 + dx;
            int py1 = y1 + dy;
            int px2 = x2 + dx;
            int py2 = y2 + dy;

            // Проверка границ изображения
            if (px1 < 0 || px1 >= width || py1 < 0 || py1 >= height ||
                px2 < 0 || px2 >= width || py2 < 0 || py2 >= height) {
                continue;
            }

            for (int c = 0; c < channels; ++c) {
                int idx1 = (py1 * width + px1) * channels + c;
                int idx2 = (py2 * width + px2) * channels + c;
                float diff = img[idx1] - img[idx2];
                sum += diff * diff;
            }
        }
    }
    return sum;
}

//Функция для создания окна вокруг пикселя
vector<unsigned char> GetWindow(const unsigned char* img, int i, int j, int width, int height, int channels, int channel, int size) {
    vector<unsigned char> window(size * size);
    for (int p = 0; p < size; p++) {
        for (int q = 0; q < size; q++) {
            int index_w = p * size + q;
            //Проверка на горизонтальные границы
            if ((i + p - size / 2 >= height) || (i + p - size / 2 < 0)) {
                //Проверка на углы
                if ((j + q - size / 2 >= width) || (j + q - size / 2 < 0)) {
                    window[index_w] = img[((i - p + size / 2) * width + j) * channels - channels * (q - size / 2) + channel];
                }
                else {
                    window[index_w] = img[((i - p + size / 2) * width + j) * channels + channels * (q - size / 2) + channel];
                }
            }
            //Проверка на вертикальные границы
            else if ((j + q - size / 2 >= width) || (j + q - size / 2 < 0)) {
                window[index_w] = img[((i + p - size / 2) * width + j) * channels - channels * (q - size / 2) + channel];
            }
            else {
                window[index_w] = img[((i + p - size / 2) * width + j) * channels + channels * (q - size / 2) + channel];
            }
        }
    }

    return window;
}

//Составление ядра для алгоритма размытия по Гауссу
vector<double> CreateGaussianKernel(int size, double sigma) {
    vector<double> kernel(size * size);
    double sum = 0.0;
    int center = size / 2;

    for (int x = -center; x <= center; x++) {
        for (int y = -center; y <= center; y++) {
            double exponent = -(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2));
            kernel[size * (x + center) + y + center] = exp(exponent) / (2 * 3.14159265 * pow(sigma, 2));
            sum += kernel[size * (x + center) + y + center];
        }
    }
    // Нормализация
    for (int k = 0; k < pow(size, 2); k++) {
        kernel[k] /= sum;
    }

    return kernel;
}

//Замена пикселей в алгоритме гауссовского размытия
vector<unsigned char> GaussianFilterAlg(const unsigned char* img, int width, int height, int channels, int size, int sigma) {
    //Создание ядра свёртки
    vector<double> kernel = CreateGaussianKernel(size, sigma);
    //Выделение памяти под новое изображение
    vector<unsigned char> modified_img(width * height * channels);
    // Обработка каждого пикселя
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Вычисляем позицию в массиве
            int index = (i * width + j) * channels;

            // Получаем значения каналов
            unsigned char r = img[index];
            unsigned char g = img[index + 1];
            unsigned char b = img[index + 2];

            vector<unsigned char> windowR = GetWindow(img, i, j, width, height, channels, 0, size);
            vector<unsigned char> windowG = GetWindow(img, i, j, width, height, channels, 1, size);
            vector<unsigned char> windowB = GetWindow(img, i, j, width, height, channels, 2, size);

            modified_img[index] = GaussPixel(windowR, kernel);  // Красный канал
            modified_img[index + 1] = GaussPixel(windowG, kernel);  // Зеленый канал
            modified_img[index + 2] = GaussPixel(windowB, kernel);  // Синий канал

            if (channels == 4) {
                modified_img[index + 3] = img[index + 3];
            }
        }
    }
    return modified_img;
}


//Замена пикселей в алгоритме фильтра по среднему значению
vector<unsigned char> AverageFilterAlg(const unsigned char* img, int width, int height, int channels) {
    vector<unsigned char> modified_img(width * height * channels);
    // Обработка каждого пикселя
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Вычисляем позицию в массиве
            int index = (i * width + j) * channels;

            // Получаем значения каналов
            unsigned char r = img[index];
            unsigned char g = img[index + 1];
            unsigned char b = img[index + 2];

            vector<unsigned char> windowR = GetWindow(img, i, j, width, height, channels, 0, 3);
            vector<unsigned char> windowG = GetWindow(img, i, j, width, height, channels, 1, 3);
            vector<unsigned char> windowB = GetWindow(img, i, j, width, height, channels, 2, 3);

            modified_img[index] = AveragePixel(windowR);  // Красный канал
            modified_img[index + 1] = AveragePixel(windowG);  // Зеленый канал
            modified_img[index + 2] = AveragePixel(windowB);  // Синий канал

            if (channels == 4) {
                modified_img[index + 3] = img[index + 3];
            }
        }
    }
    return modified_img;
}

//Замена пикселей в алгоритме медианного фильтра
vector<unsigned char> MedianFilterAlg(const unsigned char* img, int width, int height, int channels) { 
    vector<unsigned char> modified_img(width * height * channels);
    // Обработка каждого пикселя
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Вычисляем позицию в массиве
            int index = (i * width + j) * channels;

            // Получаем значения каналов
            unsigned char r = img[index];
            unsigned char g = img[index + 1];
            unsigned char b = img[index + 2];

            vector<unsigned char> windowR = GetWindow(img, i, j, width, height, channels, 0, 3);
            vector<unsigned char> windowG = GetWindow(img, i, j, width, height, channels, 1, 3);
            vector<unsigned char> windowB = GetWindow(img, i, j, width, height, channels, 2, 3);

            modified_img[index] = MedianPixel(windowR);  // Красный канал
            modified_img[index + 1] = MedianPixel(windowG);  // Зеленый канал
            modified_img[index + 2] = MedianPixel(windowB);  // Синий канал

            if (channels == 4) {
                modified_img[index + 3] = img[index + 3];
            }
        }
    }
    return modified_img;
}

// Нелокальный алгоритм сглаживания NLM
vector<unsigned char> NonLocalManualDenoise(const unsigned char* input, int width, int height, int channels) {
    vector<unsigned char> modified_img(width * height * channels);
    const float H2 = H * H;
    const int patchSize = 2 * PATCH_RADIUS + 1;
    const float normalization = 1.0f / (patchSize * patchSize * channels);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float sumWeights = 0.0f;
            vector<float> sumPixels(channels, 0.0f);

            // Поиск в окрестности
            for (int dy = -SEARCH_RADIUS; dy <= SEARCH_RADIUS; ++dy) {
                for (int dx = -SEARCH_RADIUS; dx <= SEARCH_RADIUS; ++dx) {
                    int nx = x + dx;
                    int ny = y + dy;

                    if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
                        continue;
                    }

                    // Вычисление веса
                    float dist = PatchDistance(input, width, height, channels, x, y, nx, ny);
                    float weight = exp(-dist * normalization / H2);

                    sumWeights += weight;
                    for (int c = 0; c < channels; ++c) {
                        int idx = (ny * width + nx) * channels + c;
                        sumPixels[c] += weight * input[idx];
                    }
                }
            }
            // Нормализация
            for (int c = 0; c < channels; ++c) {
                int idx = (y * width + x) * channels + c;
                modified_img[idx] = sumPixels[c] / sumWeights;
            }
        }
    }
    return modified_img;
}

//Замена пикселей в аппроксимирующем алгоритме
vector<unsigned char> ApproxiP2EPEN(const unsigned char* img, int width, int height, int channels) {
    vector<unsigned char> modified_img(width * height * channels);
    // Обработка каждого пикселя
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Вычисляем позицию в массиве
            int index = (i * width + j) * channels;

            // Получаем значения каналов
            unsigned char r = img[index];
            unsigned char g = img[index + 1];
            unsigned char b = img[index + 2];
           
            vector<unsigned char> windowR = GetWindow(img, i, j, width, height, channels, 0, 5);
            vector<unsigned char> windowG = GetWindow(img, i, j, width, height, channels, 1, 5);
            vector<unsigned char> windowB = GetWindow(img, i, j, width, height, channels, 2, 5);

            // Пример обработки - инверсия цветов
            modified_img[index] = ChangeCentralPixel(windowR);  // Красный канал
            modified_img[index + 1] = ChangeCentralPixel(windowG);  // Зеленый канал
            modified_img[index + 2] = ChangeCentralPixel(windowB);  // Синий канал

            // Если есть альфа-канал (RGBA), копируем его без изменений
            if (channels == 4) {
                modified_img[index + 3] = img[index + 3];
            }
        }
    }
    return modified_img;
}

int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "RU");
    
    // Проверка аргументов командной строки
    if (argc != 4) {
        std::cerr << "Использование: " << argv[0] << " <входной_файл> <выходной_файл> <номер_алгоритма>" << std::endl;
        std::cerr << "Номера алгоритмов:" << std::endl;
        std::cerr << "1 - Средний фильтр" << std::endl;
        std::cerr << "2 - Медианный фильтр" << std::endl;
        std::cerr << "3 - Гауссовский фильтр" << std::endl;
        std::cerr << "4 - Нелокальный алгоритм сглаживания (NLM)" << std::endl;
        std::cerr << "5 - Аппроксимация" << std::endl;
        return 1;
    }

    // Получение параметров
    const char* input_file = argv[1];
    const char* output_file = argv[2];
    int choice = std::stoi(argv[3]);

    // Загрузка изображения
    int width, height, channels;
    unsigned char* img = stbi_load(input_file, &width, &height, &channels, 0);

    if (!img) {
        std::cerr << "Ошибка загрузки изображения!" << std::endl;
        return 1;
    }

    std::cout << "Размер: " << width << "x" << height << " | Каналов: " << channels << std::endl;

    bool success = false;
    switch (choice) {
        case 1: {   // Фильтрация по среднему
            std::vector<unsigned char> modified_img = AverageFilterAlg(img, width, height, channels);
            success = stbi_write_png(output_file, width, height, channels, modified_img.data(), width * channels);
            break;
        }
        case 2: {   // Медианная фильтрация
            std::vector<unsigned char> modified_img = MedianFilterAlg(img, width, height, channels);
            success = stbi_write_png(output_file, width, height, channels, modified_img.data(), width * channels);
            break;
        }
        case 3: {   // Размытие по Гауссу
            std::vector<unsigned char> modified_img = GaussianFilterAlg(img, width, height, channels, 3, 1.5);
            success = stbi_write_png(output_file, width, height, channels, modified_img.data(), width * channels);
            break;
        }
        case 4: {   // Нелокальный алгоритм сглаживания
            std::vector<unsigned char> modified_img = NonLocalManualDenoise(img, width, height, channels);
            success = stbi_write_jpg(output_file, width, height, channels, modified_img.data(), width * channels);
            break;
        }
        case 5: {   // Алгоритм ApproxiP2EPEN
            std::vector<unsigned char> modified_img = ApproxiP2EPEN(img, width, height, channels);
            success = stbi_write_png(output_file, width, height, channels, modified_img.data(), width * channels);
            break;
        }
        default: {
            std::cerr << "Неверный номер алгоритма!" << std::endl;
            stbi_image_free(img);
            return 1;
        }
    }

    // Освобождение памяти
    stbi_image_free(img);

    if (!success) {
        std::cerr << "Ошибка сохранения изображения!" << std::endl;
        return 1;
    }

    std::cout << "Изображение успешно обработано и сохранено!" << std::endl;
    return 0;
}