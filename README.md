# ChemAssist v2

AI-powered computational chemistry toolkit for DFT, Molecular Dynamics, and more.

## 🚀 Quick Start

1. **Install dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

2. **Set up API key (for AI features):**
   
   **Option A: Environment variable (recommended):**
   ```bash
   export GROQ_API_KEY="your_api_key_here"
   ```
   
   **Option B: Streamlit secrets (for Streamlit Cloud):**
   Create `.streamlit/secrets.toml`:
   ```toml
   GROQ_API_KEY = "your_api_key_here"
   ```

3. **Run the application:**
   ```bash
   streamlit run app.py
   ```

## 🔑 API Key Setup

ChemAssist uses AI to help fix computational chemistry errors and suggest methods. You need an API key for these features:

### Quick Setup (Recommended)
1. Get a free API key from [Groq Console](https://console.groq.com/)
2. Set environment variable: `export GROQ_API_KEY="your_api_key_here"`
3. Run: `streamlit run app.py`

### Option 1: Groq (Recommended - Free Tier)
1. Visit [https://console.groq.com/](https://console.groq.com/)
2. Sign up for a free account
3. Get your API key from the dashboard
4. Set as environment variable: `export GROQ_API_KEY="your_key"`

### Option 2: OpenAI
1. Visit [https://platform.openai.com/api-keys](https://platform.openai.com/api-keys)
2. Create an account and get an API key
3. Set as environment variable: `export OPENAI_API_KEY="your_key"`

### Option 3: Streamlit Cloud Deployment
1. Deploy to Streamlit Cloud
2. Add secrets in the Streamlit Cloud dashboard:
   - Go to your app settings
   - Add `GROQ_API_KEY` or `OPENAI_API_KEY` in the secrets section

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
- Environment variables are preferred over .env files for better deployment

## 🔧 Troubleshooting

If you see "No LLM API key found" error:
1. Make sure you have a valid API key
2. Check that the environment variable is set: `echo $GROQ_API_KEY`
3. For Streamlit Cloud, verify secrets are configured in the dashboard
4. Restart the Streamlit app after setting the key 