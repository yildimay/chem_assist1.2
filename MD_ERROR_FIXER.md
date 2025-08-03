# 🚑 MD Error Fixer

## Overview
The Molecular Dynamics Error Fixer diagnoses and fixes common GROMACS simulation issues using AI. It analyzes .mdp files and .log files to identify problems and generate corrected input files.

## ✅ What's Implemented

### Core Features
- **MDP File Fixing**: Analyzes and corrects .mdp parameter files
- **Log Analysis**: Analyzes .log files for error diagnosis
- **Multi-stage Support**: Works with EM, NVT, NPT, and Production MD
- **AI-powered Diagnosis**: Uses LLM to understand and fix complex issues

### Supported Issues
- **Parameter Errors**: Invalid temperatures, pressures, cutoffs
- **File Format Issues**: Missing parameters, syntax errors
- **Performance Problems**: Memory issues, domain decomposition
- **System Issues**: Box size problems, atom overlap

### User Interface
- **File Upload**: Drag-and-drop .mdp and .log files
- **Stage Selection**: Choose simulation stage (em, nvt, npt, md)
- **Analysis Options**: Fix MDP file or analyze log only
- **Results Display**: Diagnosis, suggestions, and corrected files
- **Download Options**: Fixed .mdp and comparison files

## 🎯 How to Use

### Fix MDP File
1. **Upload files**: .mdp file and corresponding .log file
2. **Select stage**: Choose the simulation stage that failed
3. **Click "Fix MDP File"**: AI analyzes and generates corrected file
4. **Review results**: Check diagnosis and download fixed file

### Analyze Log Only
1. **Upload log file**: Just the .log file with error messages
2. **Click "Analyze Log Only"**: AI provides analysis and suggestions
3. **Review analysis**: Get detailed explanation of what went wrong

## 📁 File Support

### Input Files
- **MDP files**: .mdp, .txt
- **Log files**: .log, .out, .txt

### Output Files
- **Fixed MDP**: Corrected .mdp file ready to use
- **Comparison**: Original vs fixed file comparison
- **Analysis**: Detailed diagnosis and suggestions

## 🔧 Technical Details

### AI Analysis Process
1. **Error Detection**: Analyzes log file for error patterns
2. **Parameter Validation**: Checks .mdp file for invalid values
3. **Context Understanding**: Considers simulation stage and parameters
4. **Fix Generation**: Creates corrected .mdp file with proper parameters

### Supported Stages
- **EM**: Energy minimization
- **NVT**: Constant NVT equilibration
- **NPT**: Constant NPT equilibration
- **MD**: Production molecular dynamics

### Common Fixes
- **Temperature**: Correct negative or invalid temperature values
- **Pressure**: Fix invalid pressure coupling parameters
- **Cutoffs**: Adjust neighbor list and interaction cutoffs
- **Constraints**: Fix bond constraint settings
- **Output**: Correct output frequency parameters

## 🧪 Testing

### Test Cases
- ✅ **Parameter Errors**: Negative temperatures, invalid cutoffs
- ✅ **File Format**: Missing parameters, syntax issues
- ✅ **Performance**: Memory allocation problems
- ✅ **System Issues**: Box size and atom overlap

### Sample Issues Fixed
1. **Negative Temperature**: `ref_t = -300.0` → `ref_t = 300.0`
2. **Invalid Cutoffs**: `rlist = -1.0` → `rlist = 1.0`
3. **Missing Parameters**: Add required parameters
4. **Format Issues**: Fix syntax and formatting

## 🚀 Benefits

1. **Time Saving**: Automatic diagnosis and fixing
2. **Expert Knowledge**: AI understands GROMACS intricacies
3. **Comprehensive**: Handles multiple types of errors
4. **Educational**: Provides detailed explanations
5. **Reliable**: Generates valid, corrected files

## 📊 User Experience

### Workflow
1. **Upload Files** → **Select Options** → **Analyze** → **Download Results**

### Error Handling
- **Graceful Fallbacks**: Individual file downloads if ZIP fails
- **Clear Messages**: Helpful error messages and troubleshooting tips
- **Retry Options**: Multiple analysis types and formats

### Results Display
- **Visual Feedback**: Success/error indicators
- **File Previews**: Syntax-highlighted code display
- **Download Options**: Multiple file formats and comparisons

## 🔧 Troubleshooting

### Common Issues
1. **"LLM did not return a FIXED_MDP block"**
   - Try "Analyze Log Only" first
   - Check file formats are valid
   - Ensure log contains error messages

2. **"Invalid API Key"**
   - Run with: `./launch.sh YOUR_API_KEY`
   - Get free key from: https://console.groq.com/

3. **"No response from LLM"**
   - Check internet connection
   - Try again in a few moments

The MD Error Fixer is now fully functional and ready to help users resolve GROMACS simulation issues! 🎉 