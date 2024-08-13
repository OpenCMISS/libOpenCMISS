#include <iostream>
#include <wchar.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "opencmiss_cellml_model_f.h"
#include "opencmiss_cellml_model.hpp"

#ifdef _MSC_VER
#	include <direct.h>
#	define getcwd _getcwd
	#include <io.h>
#else
	#include <unistd.h>
#endif

int CellMLModel_CreateF(
			const char* URI,
			OpaqueCellMLObject* model
			)
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  char* inputURI = (char*)NULL;
  
  if(model && !(*model))
    {
      OpaqueCellMLObject cellml_model_definition = (OpaqueCellMLObject)NULL;
      if(URI && strlen(URI) > 1)
	{
	  // make sure we have an absolute URI
	  if(!(errorCode = getAbsoluteURI(URI, &inputURI)))
	    {
	      /*
	       * Create the model defintion object
	       */
	      CellMLModel* newModel = new CellMLModel(inputURI);
	      if(newModel)
		{
		  *model = (OpaqueCellMLObject)newModel;
		}
	      else
		{
		  errorCode = CELLML_MODEL_COULD_NOT_CREATE_MODEL_ERROR;
		}
	    }
	  if(inputURI) free(inputURI);
	}
      else
	{
	  errorCode = CELLML_MODEL_INVALID_URI_ERROR;
	}
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_NOT_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_DestroyF(
			  OpaqueCellMLObject* _ptr
			 )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* modelObj = (CellMLModel*)NULL;
  if (_ptr && *_ptr && (modelObj = (CellMLModel*)(*_ptr)))
    {
      delete modelObj;
      // want to make sure that we don't try to delete it again...
      *_ptr = (OpaqueCellMLObject)NULL;
      _ptr = (OpaqueCellMLObject*)NULL;
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }

  return(errorCode);
}

int CellMLModel_GetInitialValueF(
				 OpaqueCellMLObject _ptr,
				 const char* name,
				 double* value
				 )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getInitialValue(name,value);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetInitialValueByIndexF(
					OpaqueCellMLObject _ptr,
					const int variableType,
					const int variableIndex,
					double* value
					)
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      // need to map from 0-indexed C arrays to 1-indexed Fortran arrays
      int i = variableIndex-1;
      errorCode = model->getInitialValueByIndex(variableType,i,value);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetVariableTypeF(
				 OpaqueCellMLObject _ptr,
				 const char* name,
				 int* variableType
				 )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getVariableType(name,variableType);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetVariableIndexF(
				  OpaqueCellMLObject _ptr,
				  const char* name,
				  int* variableIndex
				  )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getVariableIndex(name,variableIndex);
      // need to map from 0-indexed C arrays to 1-indexed Fortran arrays
      if(errorCode == CELLML_MODEL_NO_ERROR) (*variableIndex) += 1;
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_SetVariableAsKnownF(
				    OpaqueCellMLObject _ptr,
				    const char* name
				    )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->setVariableAsKnown(name);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_SetVariableAsWantedF(
				     OpaqueCellMLObject _ptr,
				     const char* name
				     )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->setVariableAsWanted(name);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_InstantiateF(
			     OpaqueCellMLObject _ptr
			     )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->instantiate();
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}


int CellMLModel_GetNumberOfAlgebraicsF(
				       OpaqueCellMLObject _ptr,
				       int* numOfAlgebraics
				       )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getNumberOfAlgebraics(numOfAlgebraics);
    }
  else    
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}


int CellMLModel_GetNumberOfComputedConstantsF(
					      OpaqueCellMLObject _ptr,
					      int* numOfComputedConstants
					      )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getNumberOfComputedConstants(numOfComputedConstants);
    }
  else    
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetNumberOfConstantsF(
				      OpaqueCellMLObject _ptr,
				      int* numOfConstants
				      )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getNumberOfConstants(numOfConstants);
    }
  else    
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetNumberOfStatesF(
				   OpaqueCellMLObject _ptr,
				   int* numOfStates
				   )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getNumberOfStates(numOfStates);
    }
  else    
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int CellMLModel_GetNumberOfRatesF(
				  OpaqueCellMLObject _ptr,
				  int* numOfRates
				  )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getNumberOfStates(numOfRates);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

void CellMLModel_CallComputedConstantsRoutineF(
					       OpaqueCellMLObject _ptr,
					       double* variables
					       )
{
  //Think about errors?
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      model->callModelComputedConstantsFunction(variables);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
}
 
void CellMLModel_CallRatesRoutineF(
				   OpaqueCellMLObject _ptr,
				   double voi,
				   double* states,
				   double* rates,
				   double* wanted,
				   double* known
				   )
{
  //Think about errors?
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      model->callModelRatesFunction(voi,states,rates,wanted,known);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
}

void CellMLModel_CallVariablesRoutineF(
				       OpaqueCellMLObject _ptr,
				       double voi,
				       double* states,
				       double* rates,
				       double* variables
				       )
{
  //Think about errors?
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      model->callModelVariablesFunction(voi,states,rates,variables);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
}
 
int CellMLModel_GetLastErrorF(
			       OpaqueCellMLObject _ptr,
			       const int maxErrorStringLength,
			       char* lastError
			       )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  CellMLModel* model = (CellMLModel*)NULL;
  if (_ptr && (model = (CellMLModel*)_ptr))
    {
      errorCode = model->getLastError(maxErrorStringLength,lastError);
    }
  else
    {
      errorCode = CELLML_MODEL_PTR_IS_NULL_ERROR;
    }
  
  return(errorCode);
}

int getAbsoluteURI(
		   const char* URI,
		   char** absoluteURI
		   )
{
  int errorCode = CELLML_MODEL_NO_ERROR;
  
  if(URI)
    {
      if(absoluteURI && !*absoluteURI)
	{
	  if(strstr(URI,"://") != NULL)
	    {
	      *absoluteURI = (char*)malloc(strlen(URI)+1);
	      strcpy(*absoluteURI,URI);	  
	    }
	  else if(URI[1]==':')
	    {
	      // Windows drive letter - c:/path/to/file
	      *absoluteURI = (char*)malloc(strlen(URI)+1+7);
	      sprintf(*absoluteURI,"%s",URI);
	    }
	  else if(URI[0]=='/')
	    {
	      *absoluteURI = (char*)malloc(strlen(URI)+1+7);
	      sprintf(*absoluteURI,"file://%s",URI);
	    }
	  else
	    {
	      /* relative filename ? append absoulte path */
	      int size;
#ifdef WIN32
	      size = _MAX_PATH;
#else
	      size = pathconf(".",_PC_PATH_MAX);
#endif
	      char* cwd = (char*)malloc(size);
	      if(getcwd(cwd,size))
		{
		  *absoluteURI = (char*)malloc(strlen(cwd)+strlen(URI)+1+8);
		  sprintf(*absoluteURI,"%s/%s",cwd,URI);
		  // and make sure \'s become /'s
		  unsigned int i;
		  for(i=0; i < strlen(*absoluteURI); ++i)
		    {
		      if((*absoluteURI)[i] == '\\') (*absoluteURI)[i] = '/';
		    }
		}
	      else
		{
		  errorCode = errno;
		}
	      if(cwd) free(cwd);
	    }
	}      
    }
  
  return(errorCode);
}
