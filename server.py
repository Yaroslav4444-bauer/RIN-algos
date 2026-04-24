from flask import Flask, request, send_file
from flask_cors import CORS
import subprocess
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)
CORS(app)

# Получаем абсолютный путь к текущей директории
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Определяем пути к папкам и программе
UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
RESULT_FOLDER = os.path.join(BASE_DIR, 'results')
CPP_PROGRAM = os.path.join(BASE_DIR, 'approxymationDenoiseALG2025.exe')

# Создаем папки, если их нет
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)

# Добавляем endpoint для проверки здоровья сервера
@app.route('/health')
def health_check():
    return {'status': 'ok'}, 200

# Добавляем корневой endpoint
@app.route('/')
def home():
    return {'message': 'Server is running'}, 200

@app.route('/process', methods=['POST'])
def process():
    print("Received request")
    try:
        if 'image' not in request.files:
            print("No image in request")
            return 'No image uploaded', 400
        
        file = request.files['image']
        algo = request.form.get('algorithm', '1')
        
        if file.filename == '':
            print("No selected file")
            return 'No selected file', 400
        
        print(f"Processing file: {file.filename} with algorithm: {algo}")
        
        # Сохраняем входное изображение
        filename = secure_filename(file.filename)
        input_path = os.path.join(UPLOAD_FOLDER, filename)
        output_path = os.path.join(RESULT_FOLDER, 'result_' + filename)
        
        file.save(input_path)
        print(f"File saved to: {input_path}")
        
        # Проверяем существование C++ программы
        if not os.path.exists(CPP_PROGRAM):
            print(f"C++ program not found at: {CPP_PROGRAM}")
            raise Exception(f"C++ program not found at: {CPP_PROGRAM}")
        
        print(f"Running C++ program: {CPP_PROGRAM}")
        print(f"Input path: {input_path}")
        print(f"Output path: {output_path}")
        print(f"Algorithm: {algo}")
        
        # Запускаем C++ программу
        try:
            result = subprocess.run(
                [CPP_PROGRAM, input_path, output_path, algo],
                capture_output=True,
                text=True,
                check=True
            )
            print(f"C++ program output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            print(f"C++ program error: {e.stderr}")
            raise Exception(f"C++ program failed: {e.stderr}")
        
        # Проверяем, что выходной файл создан
        if not os.path.exists(output_path):
            raise Exception(f"Output file was not created at: {output_path}")
        
        print(f"Processing complete, result saved to: {output_path}")
        return send_file(output_path, mimetype='image/png')
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return str(e), 500

if __name__ == '__main__':
    print("Server starting on http://localhost:5000")
    print(f"Base directory: {BASE_DIR}")
    print(f"C++ program path: {CPP_PROGRAM}")
    print(f"Upload folder: {UPLOAD_FOLDER}")
    print(f"Result folder: {RESULT_FOLDER}")
    app.run(debug=True, host='0.0.0.0', port=5000)