"""Stable typed ordering shared by persisted and in-memory result tables."""
import math


def value_key(value,numeric=False):
    if value is None:return None
    if isinstance(value,float) and not math.isfinite(value):return None
    text=str(value).strip()
    if not text or text.casefold() in ('nan','none','null','undefined','n/a','na','—','--'):return None
    if numeric:
        try:
            number=float(value)
            return number if math.isfinite(number) else None
        except (TypeError,ValueError,OverflowError):return None
    return text.casefold()


def ordered(rows,key,ascending=True,numeric=False):
    """Missing values stay last in either direction; equal keys stay stable."""
    present=[];missing=[]
    for row in rows:
        value=value_key(key(row),numeric)
        if value is None:missing.append(row)
        else:present.append((value,row))
    present.sort(key=lambda item:item[0],reverse=not ascending)
    return [row for _,row in present]+missing


def is_numeric(values):
    nonempty=[value for value in values if value_key(value) is not None]
    return bool(nonempty) and all(value_key(value,True) is not None for value in nonempty)
