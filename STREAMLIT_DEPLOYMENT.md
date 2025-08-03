# Streamlit Cloud Deployment Guide

## Overview

This guide explains how to deploy ChemAssist to Streamlit Cloud using environment variables instead of `.env` files for better security and easier deployment.

## 🚀 Quick Deployment Steps

### 1. Prepare Your Repository
- ✅ Remove `.env` files (already done)
- ✅ Use environment variables (already configured)
- ✅ Update `.gitignore` (already done)

### 2. Push to GitHub
```bash
git add .
git commit -m "Configure for Streamlit Cloud deployment"
git push origin main
```

### 3. Deploy on Streamlit Cloud
1. Go to [share.streamlit.io](https://share.streamlit.io)
2. Connect your GitHub repository
3. Set deployment settings:
   - **Main file path**: `app.py`
   - **Python version**: 3.11

### 4. Configure Secrets
In your Streamlit Cloud app settings:
1. Go to "Settings" → "Secrets"
2. Add your API key:
   ```toml
   GROQ_API_KEY = "gsk_your_actual_api_key_here"
   ```

## 🔧 Environment Variable Configuration

### Local Development
```bash
# Set environment variable
export GROQ_API_KEY="your_api_key_here"

# Run the app
streamlit run app.py
```

### Streamlit Cloud
Use the secrets.toml format in the Streamlit Cloud dashboard:
```toml
GROQ_API_KEY = "gsk_your_api_key_here"
# or
OPENAI_API_KEY = "sk-your_openai_key_here"
```

## 📁 Files Changed for Deployment

### Removed Files
- ❌ `setup_api_key.py` - No longer needed
- ❌ `run_with_key.py` - No longer needed  
- ❌ `run_with_key.sh` - No longer needed
- ❌ `launch.sh` - No longer needed
- ❌ `.env` files - Replaced with environment variables

### Updated Files
- ✅ `app.py` - Removed `load_dotenv()`, updated API key messages
- ✅ `requirements.txt` - Removed `python-dotenv` dependency
- ✅ `README.md` - Updated setup instructions
- ✅ `QUICK_START.md` - Updated quick start guide
- ✅ `.gitignore` - Added Streamlit secrets and other files

### New Files
- ✅ `.streamlit/secrets.toml.example` - Example secrets configuration
- ✅ `STREAMLIT_DEPLOYMENT.md` - This deployment guide

## 🔑 API Key Setup

### Option 1: Groq (Recommended - Free)
1. Visit [console.groq.com](https://console.groq.com)
2. Sign up for free account
3. Get API key (starts with `gsk_`)
4. Add to Streamlit Cloud secrets

### Option 2: OpenAI
1. Visit [platform.openai.com](https://platform.openai.com)
2. Create account and get API key
3. Add to Streamlit Cloud secrets

## 🛠️ Testing Your Deployment

### Test API Key
```bash
# Local testing
python3 test_api_key.py YOUR_API_KEY
```

### Test App Features
1. **SMILES → 3D Modeler** - Should work without API key
2. **DFT Error Fixer** - Requires API key for AI features
3. **Input Creator** - Should work without API key

## 🚨 Troubleshooting

### "No LLM API key found"
- Check Streamlit Cloud secrets configuration
- Verify API key format (Groq: `gsk_...`, OpenAI: `sk-...`)
- Restart the app after adding secrets

### "Module not found" errors
- Check that all dependencies are in `requirements.txt`
- Verify Python version is 3.11 in Streamlit Cloud

### App not loading
- Check the main file path is set to `app.py`
- Verify repository is public or you have proper access

## 🔒 Security Best Practices

1. **Never commit API keys** to version control
2. **Use Streamlit secrets** for production deployment
3. **Use environment variables** for local development
4. **Rotate API keys** regularly
5. **Monitor API usage** to avoid unexpected charges

## 📊 Monitoring

After deployment, monitor:
- **App performance** in Streamlit Cloud dashboard
- **API usage** in your Groq/OpenAI dashboard
- **Error logs** in Streamlit Cloud logs
- **User feedback** and feature usage

## 🎉 Success!

Once deployed, your ChemAssist app will be available at:
```
https://your-app-name-your-username.streamlit.app
```

Users can access all features:
- ✅ SMILES → 3D Modeler
- ✅ DFT Input Creator & Error Fixer
- ✅ MD Input Creator & Error Fixer
- ✅ Environmental Risk Calculator 