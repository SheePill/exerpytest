from win32com.client import Dispatch

app = Dispatch("EbsOpen.Application")

print("="*80)
print("Testing: Create UniversalFluidData for EACH substance separately")
print("="*80)

# Test 1: Create UDF for substance 120 (Therminol VP1)
print("\nTest 1: Substance 120 (Therminol VP1) - DEFAULT")
ufd_120 = app.NewUniversalFluidData()
libs_120 = ufd_120.Libraries
lib_120 = libs_120.AddNewReturnObject(21)  # THERMOLIQUID

# Don't set anything - should use default (120)
h_120 = ufd_120.PropertyH_OF_PT(15.0, 320.0)
s_120 = ufd_120.PropertyS_OF_PT(15.0, 320.0)
print(f"  H = {h_120} kJ/kg")
print(f"  S = {s_120} kJ/kgK")

# Check what's active
subs_120 = lib_120.Substances
for i in range(subs_120.Count):
    sub = subs_120.Item(i)
    if sub is not None and sub.Fraction > 0:
        print(f"  Active: {sub.Substance} ({sub.SubstanceName})")

print("\n" + "="*80)
print("Test 2: Try to force substance 122 by manipulating indices")
print("="*80)

# What if we need to use GetIDsOfNames or some COM trick?
ufd_122 = app.NewUniversalFluidData()
libs_122 = ufd_122.Libraries
lib_122 = libs_122.AddNewReturnObject(21)
subs_122 = lib_122.Substances

# Try to get substance 122 directly by index (it's at index 5)
sub_122_obj = subs_122.Item(5)
print(f"Substance at index 5: {sub_122_obj.Substance} ({sub_122_obj.SubstanceName})")
print(f"  Current Fraction: {sub_122_obj.Fraction}")

# Try setting Fraction directly on the object
try:
    # First clear all
    for i in range(subs_122.Count):
        s = subs_122.Item(i)
        if s is not None:
            s.Fraction = 0.0
    
    # Now set 122 to 1.0
    sub_122_obj.Fraction = 1.0
    print(f"  After setting Fraction = 1.0: {sub_122_obj.Fraction}")
    
    # Check if it stuck
    check_sub = subs_122.Item(5)
    print(f"  Checking again: {check_sub.Fraction}")
    
    # Try calculation
    h_122 = ufd_122.PropertyH_OF_PT(15.0, 320.0)
    s_122 = ufd_122.PropertyS_OF_PT(15.0, 320.0)
    print(f"  H = {h_122} kJ/kg")
    print(f"  S = {s_122} kJ/kgK")
    
    if h_122 != h_120:
        print("  ✓ SUCCESS! Different enthalpy means substance 122 is active!")
    else:
        print("  ✗ FAILED: Same enthalpy, substance 122 not active")
        
except Exception as e:
    print(f"  Error: {e}")

print("\n" + "="*80)
print("Test 3: Check if Fraction property is read-only")
print("="*80)

ufd_test = app.NewUniversalFluidData()
libs_test = ufd_test.Libraries
lib_test = libs_test.AddNewReturnObject(21)
subs_test = lib_test.Substances

# Check if Fraction is in _prop_map_put_ (writable)
sub_test = subs_test.Item(5)
print(f"Substance object type: {type(sub_test)}")
print(f"Has Fraction setter: {'Fraction' in dir(sub_test)}")

# Try accessing the COM object properties directly
try:
    print(f"\nCOM properties (_prop_map_put_):")
    if hasattr(sub_test, '_prop_map_put_'):
        print(f"  {sub_test._prop_map_put_}")
    else:
        print("  Not accessible")
except Exception as e:
    print(f"  Error: {e}")
