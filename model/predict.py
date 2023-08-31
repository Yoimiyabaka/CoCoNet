import pandas as pd
import numpy as np
import joblib

def predict(feat: pd.DataFrame, method='ensemble', model_dir='./model', out_dir='./data'):

    vars = list(feat['Var'])

    # prepare for features
    feats = load_feats(feat)

    # standarization
    scaler = load_scalar(model_dir)
    model = load_model(method, model_dir)
    
    feats = scaler.transform(feats)
    labels, probs = predict_with_prob(feats, model, method)

    if probs is not None:
        out = pd.DataFrame({"variants": vars, "pred_labels": labels, "pred_probs": probs})
    else:
        out = pd.DataFrame({"variants": vars, "pred_labels": labels})

    out.to_csv(r'{0}/results.csv'.format(out_dir), index=False)


def predict_with_prob(data, model, method='ensemble'):
    # Check if the model has a predict_proba() method
    if hasattr(model, 'predict_proba'):
        labels = model.predict(data)
        if method == 'ensemble':
            prob_list = []
            for single_model in model.estimators_:
                if hasattr(single_model, "predict_proba"):
                    prob = single_model.predict_proba(data)[:, 1]  
                    prob_list.append(prob)
            probs = np.mean(prob_list, axis=0)
        else:
            probs = model.predict_proba(data)[:, 1]     # positive prob
        return labels.tolist(), probs.tolist()
    else:
        labels = model.predict(data)
        return labels.tolist(), None


def load_model(method, model_dir):
    if method == 'ensemble':
        model = joblib.load(r'{0}/ensemble.pkl'.format(model_dir))
    elif method == 'single':
        model = joblib.load(r'{0}/single.pkl'.format(model_dir))
    else:
        raise ValueError("Model not supported")

    return model


def load_scalar(model_dir):
    scaler = joblib.load(r'{0}/scaler.pkl'.format(model_dir))
    return scaler


def load_feats(feat:pd.DataFrame):

    feat_names = [
        'CS', 'CSS', 'P', 'D', 'B', 'C', 'Q', 'QC', 'P_Score', 
        'Sub_Cen', 'Num_cc', 'Max_cc_ratio', 'hitting_time', 'sens', 'commute_time'
    ]

    return feat[feat_names].to_numpy()