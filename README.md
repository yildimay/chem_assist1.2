# ChemAssist v2

AI-powered computational chemistry toolkit for DFT, Molecular Dynamics, and more.

## 🚀 Quick Start

1. **Install dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

2. **Set up API key (for AI features):**
   ```bash
   python3 setup_api_key.py
   ```
   Or manually edit the `.env` file:
   ```bash
   echo "GROQ_API_KEY=your_api_key_here" > .env
   ```

3. **Run the application:**
   
   **Option A: Easy launcher (recommended):**
   ```bash
   # Without API key (basic features only)
   ./launch.sh
   
   # With API key (full AI features)
   ./launch.sh YOUR_API_KEY
   ```
   
   **Option B: Python script:**
   ```bash
   python3 run_with_key.py YOUR_API_KEY
   ```
   
   **Option C: Manual run:**
   ```bash
   export PATH="/Users/yildiray/Library/Python/3.9/bin:$PATH"
   streamlit run app.py
   ```

## 🔑 API Key Setup

ChemAssist uses AI to help fix computational chemistry errors and suggest methods. You need an API key for these features:

### Quick Setup (Recommended)
1. Get a free API key from [Groq Console](https://console.groq.com/)
2. Run: `python3 run_with_key.py YOUR_API_KEY`
3. That's it! The app will start with your API key set.

### Option 1: Groq (Recommended - Free Tier)
1. Visit [https://console.groq.com/](https://console.groq.com/)
2. Sign up for a free account
3. Get your API key from the dashboard
4. Run: `python3 setup_api_key.py` and choose option 1

### Option 2: OpenAI
1. Visit [https://platform.openai.com/api-keys](https://platform.openai.com/api-keys)
2. Create an account and get an API key
3. Run: `python3 setup_api_key.py` and choose option 2

### Option 3: Manual Setup
Create a `.env` file in the project root:
```bash
GROQ_API_KEY=your_groq_api_key_here
# or
OPENAI_API_KEY=your_openai_api_key_here
```

## 🛠️ Features

- **DFT Tools:**
  - Input Creator (from PDF or SMILES)
  - Error Fixer (AI-powered diagnosis and fixes)

- **Molecular Dynamics Tools:**
  - Input Creator
  - Error Fixer

- **Other Tools:**
  - Environmental Risk Calculator
  - SMILES → 3D Modeler

## 📝 Notes

- The app will show a warning in the sidebar if no API key is configured
- Non-AI features (like input creation from templates) work without API keys
- Error fixing and AI suggestions require valid API keys

## 🔧 Troubleshooting

If you see "No LLM API key found" error:
1. Make sure you have a valid API key
2. Check that the `.env` file exists and contains your key
3. Restart the Streamlit app after adding the key 